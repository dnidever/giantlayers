"""Composition-dependent LTE continuum and line opacity on a fixed MARCS atmosphere."""
module CompositionOpacity

using Korg
using TOML

function model_columns(path)
    lines = readlines(path)
    abundance_header = findfirst(line -> startswith(line, "Logarithmic chemical number abundances"), lines)
    depth_header = findfirst(line -> occursin("Number of depth points", line), lines)
    density_header = findfirst(line -> occursin("KappaRoss", line) && occursin("Density", line), lines)
    if isnothing(abundance_header) || isnothing(depth_header) || isnothing(density_header)
        throw(ArgumentError("Expected a complete, uncompressed MARCS .mod file"))
    end
    abundances = parse.(Float64, split(join(lines[abundance_header+1:depth_header-1], " ")))
    count = parse(Int, first(split(lines[depth_header])))
    density = Float64[]
    for i in 1:count
        columns = split(lines[density_header+i])
        parse(Int, columns[1]) == i || throw(ArgumentError("MARCS density rows are out of order"))
        push!(density, parse(Float64, columns[4]))
    end
    if length(abundances) != 92 || abundances[1] != 12 || !all(isfinite, abundances)
        throw(ArgumentError("Expected 92 finite MARCS abundances with A(H) = 12"))
    end
    all(x -> isfinite(x) && x > 0, density) || throw(ArgumentError("MARCS density must be positive"))
    (; abundances, density)
end

"""Apply a metal offset and absolute per-element log epsilon overrides to a MARCS mixture."""
function composition(reference, config)
    allowed = ("metallicity_offset", "log_epsilon")
    all(key -> key in allowed, keys(config)) || throw(ArgumentError("Unknown composition setting"))
    offset = get(config, "metallicity_offset", 0.0)
    offset isa Real && isfinite(offset) || throw(ArgumentError("metallicity_offset must be finite"))

    abundances = copy(reference)
    abundances[3:end] .+= offset
    overrides = get(config, "log_epsilon", Dict())
    overrides isa AbstractDict || throw(ArgumentError("Use a [log_epsilon] table"))
    for (element, value) in overrides
        number = get(Korg.atomic_numbers, element, 0)
        number != 0 || throw(ArgumentError("Unknown element: $element"))
        value isa Real && isfinite(value) || throw(ArgumentError("Invalid abundance for $element"))
        number == 1 && value != 12 && throw(ArgumentError("Hydrogen must have A(H) = 12"))
        abundances[number] = value
    end
    return abundances
end

"""Compute opacities on a vacuum wavelength grid in Å, with microturbulence in km/s."""
function opacity(model_path, linelist, config, wavelengths; vmic=2.0)
    columns = model_columns(model_path)
    abundances = composition(columns.abundances, config)
    atmosphere = Korg.read_model_atmosphere(model_path)
    continuum = Korg.synthesize(
        atmosphere, empty(linelist), abundances, wavelengths;
        vmic, hydrogen_lines=false, return_cntm=false,
    )
    total = Korg.synthesize(
        atmosphere, linelist, abundances, wavelengths;
        vmic, hydrogen_lines=true, return_cntm=false,
        use_chemical_equilibrium_from=continuum,
    )
    # Korg returns coefficients per cm. Keep the reference atmosphere's density
    # when converting them to mass opacities for this fixed-structure calculation.
    line = max.(total.alpha .- continuum.alpha, 0) ./ columns.density
    continuum_mass = continuum.alpha ./ columns.density
    return (;
        wavelengths=total.wavelengths, line, continuum=continuum_mass,
        abundances, number_densities=total.number_densities,
        electron_density=total.electron_number_density,
    )
end

function main(args)
    if length(args) != 9
        throw(ArgumentError(
            "Usage: julia composition_opacity.jl MODEL LINELIST FORMAT COMPOSITION START_A STOP_A STEP_A VMIC OUTPUT.csv"
        ))
    end
    model, path, format, mixture, first, last, increment, microturbulence, output = args
    start, stop, step, vmic = parse.(Float64, (first, last, increment, microturbulence))
    if !all(isfinite, (start, stop, step, vmic)) || !(0 < start < stop && 0 < step <= stop-start && vmic >= 0)
        throw(ArgumentError("Use increasing positive wavelengths, a positive step, and nonnegative vmic"))
    end
    linelist = format == "apogee-dr17" ? Korg.get_APOGEE_DR17_linelist() :
               Korg.read_linelist(path; format)
    result = opacity(model, linelist, TOML.parsefile(mixture), start:step:stop; vmic)
    open(output, "w") do io
        println(io, "layer,wavelength_vacuum_A,line_opacity_cm2_g,continuum_opacity_cm2_g,electron_density_cm3")
        for (j, wavelength) in enumerate(result.wavelengths), i in axes(result.line, 1)
            println(io, join((i, wavelength, result.line[i,j], result.continuum[i,j],
                              result.electron_density[i]), ','))
        end
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CompositionOpacity.main(ARGS)
end
