"""
Creates a DataFrame with a column specifying the union of all Dict keys and one further column for each Dict.
Each row corresponds to a Dict key.
"""
function merge_dicts_to_df(; dicts...)
    out = DataFrame()
    out.key = Symbol[]

    dict_names = keys(dicts)
    dict_list = collect(values(dicts))
    dict_types = map(dict -> Union{valtype(dict), Missing}, dict_list)
    for (dict_name, dict_type) in zip(dict_names, dict_types)
        out[!, dict_name] = dict_type[]
    end

    for key in union(keys.(dict_list)...)
        push!(out, NamedTuple{(:key, dict_names...), Tuple{Symbol, dict_types...}}((key, get.(dict_list, key, missing)...)))
    end

    out
end

function apply_df_structure(df::AbstractDataFrame, structure::AbstractVector{<: Tuple{Union{Symbol, Missing}, String, Any}})
    struct_vec = map(structure) do col
        if ismissing(col[1])
            if isempty(df)
                return Symbol(col[2]) => typeof(col[3])[]
            else
                return Symbol(col[2]) => col[3]
            end
        else
            return Symbol(col[2]) => df[!, col[1]]
        end
    end

    DataFrame(struct_vec...)
end

function read_csv(path::AbstractString)
    DataFrame!(CSV.File(path, normalizenames = true))
    # names!(df, DataFrames.identifier.(string.(names(df))))
end

function read_csv(path::AbstractString, columns::AbstractVector{<: Pair{Symbol, Symbol}})
    df = read_csv(path)
    rename!(df, columns)
    df[!, last.(columns)]
end

function write_csv(path::AbstractString, df::AbstractDataFrame)
    CSV.write(path, df, escapechar = '"')
end
function write_tsv(path::AbstractString, df::AbstractDataFrame)
    CSV.write(path, df, delim = '\t', escapechar = '"')
end

function read_id_all(path; single_file::Bool = true)
    columns = [
        ("File Name", :filename, String),
        ("Protein Name", :protein, String),
        ("Peptide Sequence", :seq, String),
        ("Peptide Modified Sequence Monoisotopic Masses", :seq_semi_mod_sky, String),
        ("Modified Sequence Monoisotopic Masses", :seq_mod_sky, String),
        ("Isotope Label Type", :is_heavy, String),
        ("Precursor Mz", :pre_mz, Float64),
        ("Precursor Charge", :charge, Int),
        ("Precursor Ion Formula", :pre_formula, String),
        ("Best Retention Time", :rt, Union{Missing, Float64}),
        ("Min Start Time", :rt_min, Union{Missing, Float64}),
        ("Max End Time", :rt_max, Union{Missing, Float64}),
        ("Predicted Retention Time", :rt_predicted, Union{Missing, Float64}),
        ("Total Area MS1", :area_ms1, Union{Missing, Float64}),
        ("Total Area Fragment", :area_ms2, Union{Missing, Float64}),
        ("Isotope Dot Product", :idotp, Union{Missing, Float64}),
        ("Library Dot Product", :dotp, Union{Missing, Float64}),
        ("Ratio Dot Product", :rdotp, Union{Missing, Float64}),
        ("Detection Q Value", :q_value, Union{Missing, Float64})
    ]

    df = DataFrame!(CSV.File(path, types = Dict(n => t for (n, _, t) in columns), missingstrings = ["", "#N/A"]))

    colnames = Set(names(df))
    columns_found = filter(column -> column[1] in colnames, columns)
    rename!(df, [Symbol(pre) => post for (pre, post, _) in columns_found])

    missing_cols = setdiff(getindex.(columns, 2), propertynames(df))
    if !isempty(missing_cols)
        @warn "Not all columns of canonical id_all present: " * join(missing_cols, ", ")
    end

    df.is_heavy = df.is_heavy .== "heavy"
    # allowmissing!(df, setdiff([:rt, :rt_min, :rt_max, :area, :q_value], missing_cols))

    df
end

function parse_peptide_list(path)
    map(readlines(path)) do line
        # remove comments
        line = match(r"^(.*?)(?>#.*)?$", line).captures[1]
        if length(line) == 0
            Peptide[]
        elseif endswith(line, ".txt")
            parse_peptide_list(line)
        else
            parse_sil_sequence(line)
        end
    end |>
    x -> vcat(x...)
end

function transition_df_from_precursors(precursors::AbstractVector{PrecursorIon}; min_length = 2, max_nl_count = 3)
    transitions_df = DataFrame(pre_mod = PrecursorIon[], trans = PeptideIon[])
    for pre_mod in precursors
        transs = alltransitions(pre_mod, P_NT_y, P_CT_b; min_length = 1, precursor = true, internal = true, neutrallosses = true)
        for trans in transs
            push!(transitions_df, (pre_mod = pre_mod, trans = trans))
        end
    end
    # transitions_df = combine(:trans => only => :trans, groupby(precursors_df, Not(:trans)))

    transitions_df.iondesc = iondesc.(transitions_df.trans)
    transitions_df.fragmenttype = fragmenttype.(pepfragment.(transitions_df.trans))
    transitions_df.charge = charge.(transitions_df.trans)
    transitions_df.ordinate = ordinate.(pepfragment.(transitions_df.trans))
    transitions_df.formula = formula.(transitions_df.trans)
    transitions_df.nl_count = length.(collect.(Peptidomics.neutrallosses.(pepfragment.(transitions_df.trans))))

    filter!(row -> row.ordinate >= min_length, transitions_df)
    filter!(row -> row.nl_count <= max_nl_count, transitions_df)

    transition_prio_df = DataFrame([
        (fragmenttype = "M", nl_count = 0),
        (fragmenttype = "M", nl_count = 1),
        (fragmenttype = "y", nl_count = 0),
        (fragmenttype = "b", nl_count = 0),
        (fragmenttype = "y", nl_count = 1),
        (fragmenttype = "b", nl_count = 1),
        (fragmenttype = "yb", nl_count = 0),
        (fragmenttype = "yb", nl_count = 1),
        (fragmenttype = "M", nl_count = 2),
        (fragmenttype = "y", nl_count = 2),
        (fragmenttype = "b", nl_count = 2),
        (fragmenttype = "yb", nl_count = 2),
        (fragmenttype = "M", nl_count = 3),
        (fragmenttype = "y", nl_count = 3),
        (fragmenttype = "b", nl_count = 3),
        (fragmenttype = "yb", nl_count = 3),
    ])
    transition_prio_df.prio = 1:nrow(transition_prio_df)

    nrow_pre = nrow(transitions_df)
    transitions_df = innerjoin(transitions_df, transition_prio_df; on = [:fragmenttype, :nl_count], validate = (false, true))
    @assert nrow_pre == nrow(transitions_df)

    transitions_out_df = copy(transitions_df)
    transitions_out_df.pre_charge = charge.(transitions_out_df.pre_mod)
    transitions_out_df.pre_mz = mz.(transitions_out_df.pre_mod)
    transitions_out_df.pep_mod = peptide.(transitions_out_df.pre_mod)
    transitions_out_df.seq = string.(bare.(transitions_out_df.pep_mod))
    transitions_out_df.seq_mod = string.(clearsil.(transitions_out_df.pep_mod))
    transitions_out_df.seq_mod_sil = string.(transitions_out_df.pep_mod)

    transitions_out_df.mz = mz.(transitions_out_df.formula)
    transitions_out_df.formula = formula_to_str.(transitions_out_df.formula)

    select!(transitions_out_df, [
        :seq,
        :seq_mod,
        :seq_mod_sil,
        :pre_charge,
        :pre_mz,
        :iondesc,
        :formula,
        :fragmenttype,
        :charge,
        :ordinate,
        :nl_count,
        :prio,
        :mz
    ])

    # CSV.write(rawfile * "_transitions.csv", transitions_out_df)
    transitions_out_df
end

# allowing for a transitions and removing internals
function transition_df_from_precursors_select(precursors::AbstractVector{PrecursorIon}; ions = ['M', 'y', 'b', 'a'], min_length = 1, max_nl_count = 2)
    transitions_df = DataFrame(pre_mod = PrecursorIon[], trans = PeptideIon[])
    for pre_mod in precursors
        if 'M' in ions
            transs = alltransitions(pre_mod, P_NT_y; min_length = 1, precursor = true, neutrallosses = true)
            for trans in transs
                fragmenttype(pepfragment(trans)) != "M" && continue
                push!(transitions_df, (pre_mod = pre_mod, trans = trans))
            end
        end
        if 'y' in ions
            transs = alltransitions(pre_mod, P_NT_y; min_length = 1, precursor = false, neutrallosses = true)
            for trans in transs
                push!(transitions_df, (pre_mod = pre_mod, trans = trans))
            end
        end
        if 'b' in ions
            transs = alltransitions(pre_mod, P_CT_b; min_length = 1, precursor = false, neutrallosses = true)
            for trans in transs
                push!(transitions_df, (pre_mod = pre_mod, trans = trans))
            end
        end
        if 'a' in ions
            transs = alltransitions(pre_mod, P_CT_a; min_length = 1, precursor = false, neutrallosses = true)
            for trans in transs
                push!(transitions_df, (pre_mod = pre_mod, trans = trans))
            end
        end
    end
    # transitions_df = combine(:trans => only => :trans, groupby(precursors_df, Not(:trans)))

    transitions_df.iondesc = iondesc.(transitions_df.trans)
    transitions_df.fragmenttype = fragmenttype.(pepfragment.(transitions_df.trans))
    transitions_df.charge = charge.(transitions_df.trans)
    transitions_df.ordinate = ordinate.(pepfragment.(transitions_df.trans))
    transitions_df.formula = formula.(transitions_df.trans)
    transitions_df.nl_count = length.(collect.(Peptidomics.neutrallosses.(pepfragment.(transitions_df.trans))))

    filter!(row -> row.ordinate >= min_length, transitions_df)
    filter!(row -> row.nl_count <= max_nl_count, transitions_df)

    transition_prio_df = DataFrame([
        (fragmenttype = "M", nl_count = 0),
        (fragmenttype = "M", nl_count = 1),
        (fragmenttype = "y", nl_count = 0),
        (fragmenttype = "b", nl_count = 0),
        (fragmenttype = "a", nl_count = 0),
        (fragmenttype = "y", nl_count = 1),
        (fragmenttype = "b", nl_count = 1),
        (fragmenttype = "a", nl_count = 1),
        (fragmenttype = "M", nl_count = 2),
        (fragmenttype = "y", nl_count = 2),
        (fragmenttype = "b", nl_count = 2),
        (fragmenttype = "a", nl_count = 2),
        (fragmenttype = "M", nl_count = 3),
        (fragmenttype = "y", nl_count = 3),
        (fragmenttype = "b", nl_count = 3),
        (fragmenttype = "a", nl_count = 3),
    ])
    transition_prio_df.prio = 1:nrow(transition_prio_df)

    nrow_pre = nrow(transitions_df)
    transitions_df = innerjoin(transitions_df, transition_prio_df; on = [:fragmenttype, :nl_count], validate = (false, true))
    @assert nrow_pre == nrow(transitions_df)

    transitions_out_df = copy(transitions_df)
    transitions_out_df.pre_charge = charge.(transitions_out_df.pre_mod)
    transitions_out_df.pre_mz = mz.(transitions_out_df.pre_mod)
    transitions_out_df.pep_mod = peptide.(transitions_out_df.pre_mod)
    transitions_out_df.seq = string.(bare.(transitions_out_df.pep_mod))
    transitions_out_df.seq_mod = string.(clearsil.(transitions_out_df.pep_mod))
    transitions_out_df.seq_mod_sil = string.(transitions_out_df.pep_mod)

    transitions_out_df.mz = mz.(transitions_out_df.formula)
    transitions_out_df.formula = formula_to_str.(transitions_out_df.formula)

    select!(transitions_out_df, [
        :seq,
        :seq_mod,
        :seq_mod_sil,
        :pre_charge,
        :pre_mz,
        :iondesc,
        :formula,
        :fragmenttype,
        :charge,
        :ordinate,
        :nl_count,
        :prio,
        :mz
    ])

    # CSV.write(rawfile * "_transitions.csv", transitions_out_df)
    transitions_out_df
end
