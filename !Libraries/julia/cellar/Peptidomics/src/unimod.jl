using EzXML

import Base.convert, Base.show, Base.print, Base.+, Base.-, Base.*, Base.isequal, Base.hash

unimod_path = joinpath(dirname(pathof(Peptidomics)), "..", "data", "unimod.xml")
include_dependency(unimod_path)
xdoc = readxml(unimod_path)
# xdoc = readxml(pathof(Proteomics))

x = root(xdoc)

function showproteomicenum(io::IO, x::T) where T <: Enum
    invoke(show, Tuple{IO, Enum}, io, x)
    print(io, ": ")
    print(io, x)
end

## Elements

els = findall("/umod:unimod/umod:elements/umod:elem", x)

struct ProtElementData
    title::String
    full_name::String
    mono_mass::Float64
    avge_mass::Float64
end

export ProtElements
ProtElements = [
    ProtElementData(
        el["title"],
        el["full_name"],
        parse(Float64, el["mono_mass"]),
        parse(Float64, el["avge_mass"])
    ) for el in els
]

els_names = map(x -> x => Symbol("P_EL_" * x), getfield.(ProtElements, :title))

enum_expr = :(@enum ProtElement::UInt8)
function gen_enum_ele_entry(sym, i)
    tmplt = :(x = i)
    tmplt.args[1] = sym
    tmplt.args[2] = i
    tmplt
end
enum_entries = [:($sym = $i) for (i, (_, sym)) in enumerate(els_names)]
enum_expr.args = vcat(enum_expr.args, enum_entries)
eval(enum_expr)

export ProtElement
for (_, el) in els_names
    @eval export $el
end

export str_to_el
str_to_el = Dict(str => eval(sym) for (str, sym) in els_names)
parse(::Type{ProtElement}, x::AbstractString) = str_to_el[x]
print(io::IO, x::ProtElement) = print(io, proteomics(x).title)
show(io::IO, x::ProtElement) = showproteomicenum(io, x)
export proteomics
proteomics(el::ProtElement) = ProtElements[Integer(el)]
export proteomics_el
proteomics_el(el::AbstractString) = proteomics(parse(ProtElement, el))

export ProtFormula
ProtFormula = Dict{ProtElement, Int}
zero(::Type{ProtFormula}) = Dict{ProtElement, Int}()
+(x::Vararg{ProtFormula, N}) where N = filter(x -> x[2] != 0, merge(+, x...))
-(x::ProtFormula) = x * -1
-(x::Vararg{ProtFormula, N}) where N = filter(x -> x[2] != 0, merge(-, x...))
function *(dict::ProtFormula, x::Int)
    x == 0 && return zero(typeof(dict))
    Dict(k => v * x for (k, v) in pairs(dict))
end
function isequal(x::ProtFormula, y::ProtFormula)
    pointer_from_objref(x) == pointer_from_objref(y) && return true
    all(==(0), values(x - y))
end
function hash(x::ProtFormula, h::UInt)
    invoke(hash, Tuple{Dict, UInt}, filter(x -> x[2] != 0, x), h)
end

export mass
mass(el::ProtElement) = proteomics(el).mono_mass
mass(formula::ProtFormula) = sum(mass(k) * v for (k, v) in formula)
charge(formula::ProtFormula) = -get(formula, P_EL_e, 0)
mz(formula::ProtFormula) = mass(formula) / abs(charge(formula))

export avgmass
avgmass(el::ProtElement) = proteomics(el).avge_mass
avgmass(formula::ProtFormula) = sum(avgmass(k) * v for (k, v) in formula)

els = nothing

## AAs

aas = findall("/umod:unimod/umod:amino_acids/umod:aa", x)

struct ProtAAData
    title::String
    three_letter::String
    full_name::String
    mono_mass::Float64
    avge_mass::Float64
    formula::ProtFormula
end

function has_formula(node::EzXML.Node)
    !isnothing(findfirst("umod:element", node))
end

function parse_formula(node::EzXML.Node)
    els = findall("umod:element", node)
    [str_to_el[el["symbol"]] => parse(Int, el["number"]) for el in els]
end

export formula_to_str
formula_to_str(x::ProtFormula) = formula_to_str(collect(x))
function formula_to_str(x::AbstractVector{Pair{ProtElement, Int}})
    function elsortby(x)
        c = match(r"^(\d*)([A-Za-z]+)$", string(x)).captures
        (c[2], c[1] == "" ? 0 : parse(Int, c[1]))
    end

    reduce(*,
        # C2 -> C2, 13C13 -> [13]C13, C0 -> "", C1 -> C
        replace(string(k), r"^(\d+)" => x -> x == "" ? x : "[" * x * "]") * (v == 1 ? "" : string(v))
        # sort by element, alphabetically, main first, then isotopes
        for (k, v) in sort(x, by = x -> elsortby(x[1]))
        if v != 0
    )
end

function parse_aa(node::EzXML.Node)
    ProtAAData(
        node["title"],
        node["three_letter"],
        node["full_name"],
        parse(Float64, node["mono_mass"]),
        parse(Float64, node["avge_mass"]),
        Dict(parse_formula(node))
    )
end

export ProtAAs
ProtAAs = parse_aa.(aas)

aas_names = map(x -> x => Symbol("P_AA_" * replace(x, "-" => "_")), getfield.(ProtAAs, :title))

enum_expr = :(@enum ProtAA::UInt8)
enum_entries = [:($sym = $i) for (i, (_, sym)) in enumerate(aas_names)]
enum_expr.args = vcat(enum_expr.args, enum_entries)
eval(enum_expr)

export ProtAA
for (_, aa) in aas_names
    @eval export $aa
end

export str_to_aa
str_to_aa = Dict(x => eval(sym) for (x, sym) in aas_names)
parse(::Type{ProtAA}, x::AbstractString) = str_to_aa[x]
print(io::IO, x::ProtAA) = print(io, proteomics(x).title)
show(io::IO, x::ProtAA) = showproteomicenum(io, x)
proteomics(aa::ProtAA) = ProtAAs[Integer(aa)]
export proteomics_aa
proteomics_aa(aa::AbstractString) = proteomics(str_to_aa[aa])

export fullsil
function fullsil(aa::ProtAA)
    #Label:13C(6)15N(4)
    formula = proteomics(aa).formula
    ccount = get(formula, P_EL_C, 0)
    ncount = get(formula, P_EL_N, 0)
    if ncount == 0 && ccount == 0 error() end

    clabel = ccount > 0 ? "13C($ccount)" : ""
    nlabel = ncount > 0 ? "15N($ncount)" : ""
    parse(UnimodLabel, "Label:$clabel$nlabel")
end

function mass(aa::ProtAA) proteomics(aa).mono_mass end
function avgmass(aa::ProtAA) proteomics(aa).avge_mass end
export formula
function formula(aa::ProtAA) proteomics(aa).formula end

aas = nothing

## AANeutralLosses

struct AANeutralLoss
    title::String
    specificity::Vector{ProtAA}
    mono_mass::Float64
    avge_mass::Float64
    formula::ProtFormula
end
function newloss(specificity::AbstractVector{ProtAA}, formula_pairs::AbstractVector{Pair{ProtElement, Int}})
    formula = -Dict(formula_pairs)
    AANeutralLoss('-' * formula_to_str(formula_pairs), specificity, mass(formula), avgmass(formula), formula)
end
AANeutralLosses = [
    newloss([P_AA_S, P_AA_T, P_AA_D, P_AA_E], [P_EL_H => 2, P_EL_O => 1]),
    newloss([P_AA_N, P_AA_Q, P_AA_K, P_AA_R], [P_EL_N => 1, P_EL_H => 3]),
]

export neutraloss
neutralloss = Dict(aa => aaloss for aaloss in AANeutralLosses for aa in aaloss.specificity)
export str_to_neutraloss
str_to_neutraloss = Dict(aaloss.title => aaloss for aaloss in AANeutralLosses)
print(io::IO, x::AANeutralLoss) = print(io, x.title)

function mass(nl::AANeutralLoss) nl.mono_mass end
function avgmass(nl::AANeutralLoss) nl.avge_mass end
function formula(nl::AANeutralLoss) nl.formula end

## Unimods

mods = findall("/umod:unimod/umod:modifications/umod:mod", x)

struct UnimodNeutralLoss
    title::String
    mono_mass::Float64
    avge_mass::Float64
    formula::ProtFormula
end
@enum UnimodPosition::UInt8 Anywhere ProteinNterm AnyNterm AnyCterm ProteinCterm
# str_to_pos = Dict(string(v) => UnimodPosition(k) for (k, v) in pairs(Base.Enums.namemap(UnimodPosition)))
str_to_pos = Dict(string(pos) => pos for pos in instances(UnimodPosition))
struct UnimodSpecificity
    classification::String
    site::ProtAA
    position::UnimodPosition
    neutralloss::Union{UnimodNeutralLoss, Missing}
end
abstract type UnimodBase end
struct Unimod <: UnimodBase
    id::Int
    title::String
    full_name::String
    mono_mass::Float64
    avge_mass::Float64
    formula::ProtFormula
    specificities::Dict{Tuple{ProtAA, UnimodPosition}, UnimodSpecificity}
end
struct UnimodLabel <: UnimodBase
    id::Int
    title::String
    full_name::String
    mono_mass::Float64
    avge_mass::Float64
    formula::ProtFormula
end
struct SpecifiedUnimod
    mod::Unimod
    specificity::UnimodSpecificity
end

function parse_neutralloss(node::EzXML.Node)
    formula = parse_formula(node)
    UnimodNeutralLoss(
        '-' * formula_to_str(formula),
        -parse(Float64, node["mono_mass"]),
        -parse(Float64, node["avge_mass"]),
        -Dict(formula)
    )
end

function parse_specificity(node::EzXML.Node)
    neutrallosses = findall("umod:NeutralLoss", node)
    filter!(has_formula, neutrallosses)
    if isempty(neutrallosses)
        neutralloss = missing
    else
        neutralloss = parse_neutralloss(only(neutrallosses))
    end

    UnimodSpecificity(
        node["classification"],
        parse(ProtAA, node["site"]),
        str_to_pos[replace(node["position"], r"[\s-]" => "")],
        neutralloss
    )
end

function is_formula_sil(form::ProtFormula)
    # sil are only, where all entries are heavy and corresponding -light
    heavylight = [
        P_EL_13C => P_EL_C,
        P_EL_15N => P_EL_N,
        P_EL_18O => P_EL_O,
        P_EL_2H => P_EL_H,
    ]
    hlv = vcat(first.(heavylight), last.(heavylight))
    isempty(form) && return false
    for (h, l) in heavylight
        get(form, h, 0) == -get(form, l, 0) ||
            return false
    end
    for (k, v) in form
        k ∉ hlv && v != 0 &&
            return false
    end
    return true
end

function parse_mod(node::EzXML.Node)
    id = parse(Int32, node["record_id"])
    title = node["title"]
    full_name = node["full_name"]
    delta_nd = findfirst("umod:delta", node)
    mono_mass = parse(Float64, delta_nd["mono_mass"])
    avge_mass = parse(Float64, delta_nd["avge_mass"])
    formula = Dict(parse_formula(delta_nd))

    if startswith(title, "Label:") && is_formula_sil(formula)
        UnimodLabel(
            id,
            title,
            full_name,
            mono_mass,
            avge_mass,
            formula
        )
    else
        specificities = collect(parse_specificity.(findall("umod:specificity", node)))
        if !allunique((s.site, s.position) for s in specificities) error("Duplicate Specificity") end

        Unimod(
            id,
            title,
            full_name,
            mono_mass,
            avge_mass,
            formula,
            Dict((s.site, s.position) => s for s in specificities)
        )
    end
end

export Unimods
mods_tmp = collect(parse_mod.(mods))
Unimods = Vector{Union{UnimodBase, Missing}}(missing, maximum(x -> x.id, mods_tmp))
for mod in mods_tmp
    if !ismissing(Unimods[mod.id]) error("Duplicate Unimod ID: ", mod.id) end
    Unimods[mod.id] = mod
end
mods_tmp = nothing

# unimod_specs = Dict{Tuple{Int, ProtAA, UnimodPosition}, SpecifiedUnimod}()
# for mod in Unimods, spec in mod.specificities
#     key = (mod.id, spec.site, spec.position)
#     if haskey(unimod_specs, key) error(key) end
#     unimod_specs[key] = SpecifiedUnimod(mod, spec)
# end

# export specify
# function specify(mod::Unimod, aa::ProtAA)
#     try
#         unimod_specs[(mod.id, aa, Anywhere)]
#     catch
#         error("Unimod ", string(mod), " is not specified for ", string(aa))
#     end
# end
export specify
function specify(mod::Unimod, aa::ProtAA)
    SpecifiedUnimod(mod, Unimods[mod.id].specificities[(aa, Anywhere)])
end
export str_to_unimod
str_to_unimod = Dict(mod.title => mod for mod in Unimods if !ismissing(mod))
parse(::Type{Unimod}, x::AbstractString) = str_to_unimod[x]::Unimod
parse(::Type{UnimodLabel}, x::AbstractString) = str_to_unimod[x]::UnimodLabel
print(io::IO, x::UnimodNeutralLoss) = print(io, x.title)
function print(io::IO, x::UnimodBase; asmass = false)
    if asmass
        print(io, sign(mass(x)) >= 0 ? '+' : '-', round(Int, mass(x)))
    else
        print(io, x.title)
    end
end
print(io::IO, x::SpecifiedUnimod; asmass = false) = print(io, x.mod; asmass = asmass)
# show(io::IO, x::Unimod) = showproteomicenum(io, x)

function mass(nl::UnimodNeutralLoss) nl.mono_mass end
function mass(mod::UnimodBase) mod.mono_mass end
function mass(mod::SpecifiedUnimod) mass(mod.mod) end
function avgmass(nl::UnimodNeutralLoss) nl.avge_mass end
function avgmass(mod::UnimodBase) mod.avge_mass end
function avgmass(mod::SpecifiedUnimod) avgmass(mod.mod) end
function formula(nl::UnimodNeutralLoss) nl.formula end
function formula(mod::UnimodBase) mod.formula end
function formula(mod::SpecifiedUnimod) formula(mod.mod) end

mods = nothing

## NTerms

struct ProtNTermData
    title::String
    mono_mass::Float64
    avge_mass::Float64
    formula::Dict{ProtElement, Int}
end

ProtNTerms = Vector{ProtNTermData}()

aanterm = proteomics(P_AA_N_term)
push!(ProtNTerms, ProtNTermData(
    "H",
    aanterm.mono_mass,
    aanterm.avge_mass,
    aanterm.formula
))
# x formula
tmpform = ProtFormula(
    P_EL_C => 1,
    P_EL_O => 1,
    P_EL_H => -1
)
push!(ProtNTerms, ProtNTermData(
    "x",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))
# y formula
tmpform = ProtFormula(
    P_EL_H => 1
)
push!(ProtNTerms, ProtNTermData(
    "y",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))
# z formula
tmpform = ProtFormula(
    P_EL_N => -1,
    P_EL_H => -2
)
push!(ProtNTerms, ProtNTermData(
    "z",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))

nterms_names = map(x -> x => Symbol("P_NT_" * replace(x, "-" => "_")), getfield.(ProtNTerms, :title))

enum_expr = :(@enum ProtNTerm::UInt8)
enum_entries = [:($sym = $i) for (i, (_, sym)) in enumerate(nterms_names)]
enum_expr.args = vcat(enum_expr.args, enum_entries)
eval(enum_expr)

export ProtNTerm
for (_, nterm) in nterms_names
    @eval export $nterm
end

function specify(mod::Unimod, aa::ProtAA, ::Type{ProtNTerm})
    specs = Unimods[mod.id].specificities
    trys = [
        (aa, AnyNterm),
        (aa, ProteinNterm),
        (aa, Anywhere),
    ]
    for k in trys
        try return SpecifiedUnimod(mod, specs[k])
        catch end
    end
    error()
end
export str_to_nterm
str_to_nterm = Dict(x => eval(sym) for (x, sym) in nterms_names)
parse(::Type{ProtNTerm}, x::AbstractString) = str_to_nterm[x]
print(io::IO, x::ProtNTerm) = print(io, proteomics(x).title)
show(io::IO, x::ProtNTerm) = showproteomicenum(io, x)
function proteomics(nterm::ProtNTerm) ProtNTerms[Integer(nterm)] end
export proteomics_nterm
function proteomics_nterm(nterm::AbstractString) proteomics(str_to_nterm[nterm]) end

function mass(nterm::ProtNTerm) proteomics(nterm).mono_mass end
function avgmass(nterm::ProtNTerm) proteomics(nterm).avge_mass end
function formula(nterm::ProtNTerm) proteomics(nterm).formula end

## CTerms

struct ProtCTermData
    title::String
    mono_mass::Float64
    avge_mass::Float64
    formula::Dict{ProtElement, Int}
end

ProtCTerms = Vector{ProtCTermData}()

aacterm = proteomics(P_AA_C_term)
push!(ProtCTerms, ProtCTermData(
    "OH",
    aacterm.mono_mass,
    aacterm.avge_mass,
    aacterm.formula
))
# a formula
tmpform = ProtFormula(
    P_EL_C => -1,
    P_EL_O => -1,
    P_EL_H => -1
)
push!(ProtCTerms, ProtCTermData(
    "a",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))
# b formula
tmpform = ProtFormula(
    P_EL_H => -1
)
push!(ProtCTerms, ProtCTermData(
    "b",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))
# c formula
tmpform = ProtFormula(
    P_EL_N => 1,
    P_EL_H => 2
)
push!(ProtCTerms, ProtCTermData(
    "c",
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))

cterms_names = map(x -> x => Symbol("P_CT_" * replace(x, "-" => "_")), getfield.(ProtCTerms, :title))

enum_expr = :(@enum ProtCTerm::UInt8)
enum_entries = [:($sym = $i) for (i, (_, sym)) in enumerate(cterms_names)]
enum_expr.args = vcat(enum_expr.args, enum_entries)
eval(enum_expr)

export ProtCTerm
for (_, cterm) in cterms_names
    @eval export $cterm
end

function specify(mod::Unimod, aa::ProtAA, ::Type{ProtCTerm})
    specs = Unimods[mod.id].specificities
    trys = [
        (aa, AnyCterm),
        (aa, ProteinCterm),
        (aa, Anywhere)
    ]
    for k in trys
        try return SpecifiedUnimod(mod, specs[k])
        catch end
    end
    error()
end
export str_to_cterm
str_to_cterm = Dict(x => eval(sym) for (x, sym) in cterms_names)
parse(::Type{ProtCTerm}, x::AbstractString) = str_to_cterm[x]
print(io::IO, x::ProtCTerm) = print(io, proteomics(x).title)
show(io::IO, x::ProtCTerm) = showproteomicenum(io, x)
function proteomics(cterm::ProtCTerm) ProtCTerms[Integer(cterm)] end
export proteomics_cterm
function proteomics_cterm(cterm::AbstractString) proteomics(str_to_cterm[cterm]) end

function mass(cterm::ProtCTerm) proteomics(cterm).mono_mass end
function avgmass(cterm::ProtCTerm) proteomics(cterm).avge_mass end
function formula(cterm::ProtCTerm) proteomics(cterm).formula end

## Adducts

struct ProtAdductData
    title::String
    charge::Int
    mono_mass::Float64
    avge_mass::Float64
    formula::Dict{ProtElement, Int}
end

ProtAdducts = Vector{ProtAdductData}()

tmpform = ProtFormula(
    P_EL_H => 1,
    P_EL_e => -1
)
push!(ProtAdducts, ProtAdductData(
    "H+",
    1,
    mass(tmpform),
    avgmass(tmpform),
    tmpform
))

adduct_names = map(x -> x => Symbol("P_ADD_" * replace(replace(x, "+" => "p"), "-" => "m")), getfield.(ProtAdducts, :title))

enum_expr = :(@enum ProtAdduct::UInt8)
enum_entries = [:($sym = $i) for (i, (_, sym)) in enumerate(adduct_names)]
enum_expr.args = vcat(enum_expr.args, enum_entries)
eval(enum_expr)

export ProtAdduct
for (_, adduct) in adduct_names
    @eval export $adduct
end

export str_to_adduct
str_to_adduct = Dict(x => eval(sym) for (x, sym) in adduct_names)
parse(::Type{ProtAdduct}, x::AbstractString) = str_to_adduct[x]
print(io::IO, x::ProtAdduct) = print(io, proteomics(x).title)
show(io::IO, x::ProtAdduct) = showproteomicenum(io, x)
export proteomics_adduct
function proteomics(adduct::ProtAdduct) ProtAdducts[Integer(adduct)] end
function proteomics_adduct(adduct::AbstractString) proteomics(str_to_adduct[adduct]) end

function mass(adduct::ProtAdduct) proteomics(adduct).mono_mass end
function avgmass(adduct::ProtAdduct) proteomics(adduct).avge_mass end
function formula(adduct::ProtAdduct) proteomics(adduct).formula end

tmpform = nothing
