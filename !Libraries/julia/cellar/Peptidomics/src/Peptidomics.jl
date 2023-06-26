module Peptidomics

using BioSymbols, Combinatorics, Accessors, PrecompileTools
import Base.parse, Base.show, Base.copy, Base.isequal, Base.hash, Base.==

include("unimod.jl")

struct PeptidePos
    aa::ProtAA
    label::Union{UnimodLabel, Missing}
    mod::Union{SpecifiedUnimod, Missing}
end

export PeptideSequence
PeptideSequence = Vector{PeptidePos}

struct PeptideNTerm
    nterm::ProtNTerm
    mod::Union{SpecifiedUnimod, Missing}
end
struct PeptideCTerm
    cterm::ProtCTerm
    mod::Union{SpecifiedUnimod, Missing}
end

abstract type PeptideBase end
export Peptide
struct Peptide <: PeptideBase
    sequence::PeptideSequence
    nterm::PeptideNTerm
    cterm::PeptideCTerm
end
export PeptideFragment
abstract type PeptideFragment <: PeptideBase end
export PeptideFragmentFull
struct PeptideFragmentFull <: PeptideFragment
    peptide::Peptide
    neutrallosses::Vector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}}
    neutralloss_n::Union{UnimodNeutralLoss, Missing}
    neutralloss_c::Union{UnimodNeutralLoss, Missing}
end
export PeptideFragmentNTerm
struct PeptideFragmentNTerm <: PeptideFragment
    peptide::Peptide
    cterm::ProtCTerm
    fragment_pos::Int
    neutrallosses::Vector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}}
    neutralloss_n::Union{UnimodNeutralLoss, Missing}
end
export PeptideFragmentCTerm
struct PeptideFragmentCTerm <: PeptideFragment
    peptide::Peptide
    nterm::ProtNTerm
    fragment_pos::Int
    neutrallosses::Vector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}}
    neutralloss_c::Union{UnimodNeutralLoss, Missing}
end
export PeptideFragmentInternal
struct PeptideFragmentInternal <: PeptideFragment
    peptide::Peptide
    nterm::ProtNTerm
    cterm::ProtCTerm
    fragment_pos_n::Int
    fragment_pos_c::Int
    neutrallosses::Vector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}}
end

==(x::Peptide, y::Peptide) = isequal(x, y)
function isequal(x::Peptide, y::Peptide)
    isequal(x.sequence, y.sequence) &&
    isequal(x.nterm, y.nterm) &&
    isequal(x.cterm, y.cterm)
end
function hash(x::Peptide, h::UInt)
    reduce(xor, hash.((x.sequence, x.nterm, x.cterm), h))
end

export peptide
peptide(x::PeptideBase) = x
peptide(x::PeptideFragment) = x.peptide

export sequence
sequence(x::PeptideBase) = x.sequence
sequence(x::PeptideFragmentFull) = sequence(peptide(x))
sequence(x::PeptideFragmentNTerm) = view(sequence(peptide(x)), 1:(x.fragment_pos - 1))
sequence(x::PeptideFragmentCTerm) = view(sequence(peptide(x)), x.fragment_pos:length(sequence(peptide(x))))
sequence(x::PeptideFragmentInternal) = view(sequence(peptide(x)), x.fragment_pos_n:(x.fragment_pos_c - 1))

sequencemapper(::PeptideBase) = function(x) x end
sequencemapper(p::PeptideFragmentCTerm) = function(x) x + p.fragment_pos - 1 end
sequencemapper(p::PeptideFragmentInternal) = function(x) x + p.fragment_pos_n - 1 end

export nterm
nterm(x::PeptideBase) = x.nterm
nterm(x::PeptideFragmentFull) = nterm(peptide(x))
nterm(x::PeptideFragmentNTerm) = nterm(peptide(x))

export cterm
cterm(x::PeptideBase) = x.cterm
cterm(x::PeptideFragmentFull) = cterm(peptide(x))
cterm(x::PeptideFragmentCTerm) = cterm(peptide(x))

nlsequence(p::PeptideFragment) = p.neutrallosses

nlnterm(p::PeptideFragment) = p.neutralloss_n
nlnterm(p::PeptideFragmentCTerm) = missing
nlnterm(p::PeptideFragmentInternal) = missing

nlcterm(p::PeptideFragment) = p.neutralloss_c
nlcterm(p::PeptideFragmentNTerm) = missing
nlcterm(p::PeptideFragmentInternal) = missing

neutrallosses(p::PeptideFragment) = skipmissing(vcat(nlsequence(p), nlnterm(p), nlcterm(p)))

export PeptideIon
abstract type PeptideIon end
export PrecursorIon
struct PrecursorIon <: PeptideIon
    peptide::Peptide
    adduct::ProtAdduct
    adduct_count::Int
end
export FragmentIon
struct FragmentIon <: PeptideIon
    peptide::PeptideFragment
    adduct::ProtAdduct
    adduct_count::Int
end

peptide(x::PeptideIon) = peptide(x.peptide)
export pepfragment
pepfragment(x::PeptideIon) = x.peptide
sequence(x::PeptideIon) = sequence(x.peptide)
nterm(x::PeptideIon) = nterm(x.peptide)
cterm(x::PeptideIon) = cterm(x.peptide)
export mz
mz(pi::PeptideIon) = mass(pi) / abs(charge(pi))
export charge
charge(pi::PeptideIon) = proteomics(pi.adduct).charge * pi.adduct_count

export ionize
function ionize(p::Peptide, charge::Int, adduct::ProtAdduct=P_ADD_Hp)
    PrecursorIon(
        p,
        adduct,
        charge
    )
end
function ionize(p::PeptideFragment, charge::Int, adduct::ProtAdduct=P_ADD_Hp)
    FragmentIon(
        p,
        adduct,
        charge
    )
end

boolmap(x::AbstractVector{Vector{T}}) where T = map(x -> Vector(falses(length(x))), x)

export fragment
function fragment(p::Peptide)
    PeptideFragmentFull(
        p,
        Vector{Missing}(missing, length(sequence(p))),
        missing,
        missing
    )
end
function fragment(p::Peptide, type::ProtNTerm, pos::Int)
    seqlen = length(sequence(p))
    1 <= pos <= seqlen || error()

    PeptideFragmentCTerm(
        p,
        type,
        pos,
        Vector{Missing}(missing, seqlen - pos + 1),
        missing
    )
end
function fragment(p::Peptide, type::ProtCTerm, pos::Int)
    1 < pos <= (length(sequence(p)) + 1) || error()

    PeptideFragmentNTerm(
        p,
        type,
        pos,
        Vector{Missing}(missing, pos - 1),
        missing
    )
end
function fragment(p::Peptide, type_n::ProtNTerm, type_c::ProtCTerm, pos_n::Int, pos_c::Int)
    1 <= pos_n < pos_c <= (length(sequence(p)) + 1) || error()

    PeptideFragmentInternal(
        p,
        type_n,
        type_c,
        pos_n,
        pos_c,
        Vector{Missing}(missing, pos_c - pos_n)
    )
end
function fragment(pf::PeptideFragmentFull, type::ProtNTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 <= pos <= length(sequence(pf)) || error()

    PeptideFragmentCTerm(
        p,
        type,
        sm(pos),
        nlsequence(pf)[pos:end],
        nlnterm(pf)
    )
end
function fragment(pf::PeptideFragmentNTerm, type::ProtNTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 <= pos <= length(sequence(pf)) || error()

    PeptideFragmentInternal(
        p,
        type,
        cterm(pf),
        sm(pos),
        pf.fragment_pos,
        nlsequence(pf)[pos:end],
    )
end
function fragment(pf::PeptideFragmentCTerm, type::ProtNTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 <= pos <= length(sequence(pf)) || error()

    PeptideFragmentCTerm(
        p,
        type,
        sm(pos),
        nlsequence(pf)[pos:end],
        nlcterm(pf)
    )
end
function fragment(pf::PeptideFragmentInternal, type::ProtNTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 <= pos <= length(sequence(pf)) || error()

    PeptideFragmentInternal(
        p,
        type,
        cterm(pf),
        sm(pos),
        pf.fragment_pos_c,
        nlsequence(pf)[pos:end]
    )
end
function fragment(pf::PeptideFragmentFull, type::ProtCTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 < pos <= (length(sequence(pf)) + 1) || error()

    PeptideFragmentNTerm(
        p,
        type,
        sm(pos),
        nlsequence(pf)[1:(pos - 1)],
        nlcterm(pf)
    )
end
function fragment(pf::PeptideFragmentNTerm, type::ProtCTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 < pos <= length(sequence(pf)) || error()

    PeptideFragmentNTerm(
        p,
        type,
        sm(pos),
        nlsequence(pf)[1:(pos - 1)],
        nlnterm(pf)
    )
end
function fragment(pf::PeptideFragmentCTerm, type::ProtCTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 < pos <= length(sequence(pf)) || error()

    PeptideFragmentInternal(
        p,
        nterm(pf),
        type,
        pf.fragment_pos,
        sm(pos),
        nlsequence(pf)[1:(pos - 1)],
    )
end
function fragment(pf::PeptideFragmentInternal, type::ProtCTerm, pos::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 < pos <= length(sequence(pf)) || error()

    PeptideFragmentInternal(
        p,
        nterm(pf),
        type,
        pf.fragment_pos_n,
        sm(pos),
        nlsequence(pf)[1:(pos - 1)],
    )
end
function fragment(pf::PeptideFragment, type_n::ProtNTerm, type_c::ProtCTerm, pos_n::Int, pos_c::Int)
    p = peptide(pf)
    sm = sequencemapper(pf)
    1 <= pos_n < pos_c <= length(sequence(pf)) || error()

    PeptideFragmentInternal(
        p,
        type_n,
        type_c,
        sm(pos_n),
        sm(pos_c),
        nlsequence(pf)[pos_n:(pos_c - 1)]
    )
end

function gen_neutrallosses(sequence::AbstractVector{PeptidePos}, preset_losses::AbstractVector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}} = Vector{Union{AANeutralLoss, UnimodNeutralLoss, Missing}}(missing, length(sequence)))
    ignore_positions = .!ismissing.(preset_losses)
    nl_mods = Vector{Tuple{Int, Union{AANeutralLoss, UnimodNeutralLoss}}}()

    for (i, ppos) in enumerate(sequence)
        ignore_positions[i] && continue

        mod = ppos.mod
        if !ismissing(mod)
            nl = mod.specificity.neutralloss
            ismissing(nl) || push!(nl_mods, (i, nl))
        else
            try push!(nl_mods, (i, neutralloss[ppos.aa]))
            catch end
        end
    end

    map(powerset(nl_mods)) do nl
        nl_vec = copy(preset_losses)
        for (i, mod) in nl
            nl_vec[i] = mod
        end
        nl_vec
    end
end
# function gen_seq_neutrallosses(sequence::AbstractVector{PeptidePos}) gen_seq_neutrallosses(sequence, Int32[])

export allneutrallosses
function allneutrallosses(pf::PeptideFragmentFull)
    ntermmod = nterm(pf).mod
    do_nterm = ismissing(nlnterm(pf)) && !ismissing(ntermmod) && !ismissing(ntermmod.specificity.neutralloss)
    ctermmod = cterm(pf).mod
    do_cterm = ismissing(nlcterm(pf)) && !ismissing(ctermmod) && !ismissing(ctermmod.specificity.neutralloss)

    lossfragmentmaker(nlnterm, nlcterm) = nls -> PeptideFragmentFull(peptide(pf), nls, nlnterm, nlcterm)
    lossfragmentmakers = [lossfragmentmaker(missing, missing)]
    do_nterm && push!(lossfragmentmakers, lossfragmentmaker(ntermmod.specificity.neutralloss, missing))
    do_cterm && push!(lossfragmentmakers, lossfragmentmaker(missing, ctermmod.specificity.neutralloss))
    do_nterm && do_cterm && push!(lossfragmentmakers, lossfragmentmaker(ntermmod.specificity.neutralloss, ctermmod.specificity.neutralloss))

    (lossfragmentmaker(nls) for lossfragmentmaker in lossfragmentmakers for nls in gen_neutrallosses(sequence(pf), nlsequence(pf)))
end
function allneutrallosses(pf::PeptideFragmentNTerm)
    ntermmod = nterm(pf).mod
    do_nterm = ismissing(nlnterm(pf)) && !ismissing(ntermmod) && !ismissing(ntermmod.specificity.neutralloss)

    lossfragmentmaker(nlnterm) = nls -> PeptideFragmentNTerm(peptide(pf), cterm(pf), pf.fragment_pos, nls, nlnterm)
    lossfragmentmakers = [lossfragmentmaker(missing)]
    do_nterm && push!(lossfragmentmakers, lossfragmentmaker(ntermmod.specificity.neutralloss))

    (lossfragmentmaker(nls) for lossfragmentmaker in lossfragmentmakers for nls in gen_neutrallosses(sequence(pf), nlsequence(pf)))
end
function allneutrallosses(pf::PeptideFragmentCTerm)
    ctermmod = cterm(pf).mod
    do_cterm = ismissing(nlcterm(pf)) && !ismissing(ctermmod) && !ismissing(ctermmod.specificity.neutralloss)

    lossfragmentmaker(nlcterm) = nls ->  PeptideFragmentCTerm(peptide(pf), nterm(pf), pf.fragment_pos, nls, nlcterm)
    lossfragmentmakers = [lossfragmentmaker(missing)]
    do_cterm && push!(lossfragmentmakers, lossfragmentmaker(ctermmod.specificity.neutralloss))

    (lossfragmentmaker(nls) for lossfragmentmaker in lossfragmentmakers for nls in gen_neutrallosses(sequence(pf), nlsequence(pf)))
end
function allneutrallosses(pf::PeptideFragmentInternal)
    (PeptideFragmentInternal(peptide(pf), nterm(pf), cterm(pf), pf.fragment_pos_n, pf.fragment_pos_c, nls) for nls in gen_neutrallosses(sequence(pf), nlsequence(pf)))
end

export allfragments
function allfragments(p::PeptideBase, nterm::ProtNTerm, cterm::ProtCTerm; min_length::Int = 1, precursor::Bool = false, internal::Bool = false)
    seqlen = length(sequence(p))
    0 < min_length < seqlen || error()

    frags = Vector{PeptideFragment}()
    precursor && push!(frags, p isa PeptideFragment ? p : fragment(p))

    for i in seqlen:-1:(min_length + 1)
        push!(frags, fragment(p, nterm, seqlen - i + 2))
        push!(frags, fragment(p, cterm, i))
    end
    if (internal)
        for len in (seqlen - 2):-1:min_length, pos_from_n in 2:(seqlen - len)
            push!(frags, fragment(p, nterm, cterm, pos_from_n, pos_from_n + len))
        end

        # seqlen_n_max = seqlen - 1
        # seqlen_c_max = seqlen - 1
        # seqlen_n_min = 1
        # seqlen_c_min = 1
        # # # e.g. y ions are identical to normal n term so first internal ion is exactly like common b ion and is obsolete
        # # if formula(nterm) == formula(P_AA_N_term)
        # #     seqlen_n_min += 1
        # # end
        # # if formula(cterm) == formula(P_AA_C_term)
        # #     seqlen_c_max -= 1
        # # end
        # for i_from_n in seqlen_n_min:seqlen_n_max, i_to_c in seqlen_c_max:-1:seqlen_c_min
        #     if i_to_c - i_from_n + 1 >= min_length
        #         push!(frags, fragment(p, nterm, cterm, i_from_n, i_to_c))
        #     end
        # end
    end
    frags
end
function allfragments(p::PeptideBase, nterm::ProtNTerm; min_length::Int = 1, precursor::Bool = false)
    seqlen = length(sequence(p))
    0 < min_length < seqlen || error()

    frags = Vector{PeptideFragment}()
    precursor && push!(frags, p isa PeptideFragment ? p : fragment(p))

    for i in seqlen:-1:(min_length + 1)
        push!(frags, fragment(p, nterm, seqlen - i + 2))
    end
    frags
end
function allfragments(p::PeptideBase, cterm::ProtCTerm; min_length::Int = 1, precursor::Bool = false)
    seqlen = length(sequence(p))
    0 < min_length < seqlen || error()

    frags = Vector{PeptideFragment}()
    precursor && push!(frags, p isa PeptideFragment ? p : fragment(p))

    for i in seqlen:-1:(min_length + 1)
        push!(frags, fragment(p, cterm, i))
    end
    frags
end

function allions_helper(frags::AbstractVector{PeptideFragment}, c_max::Int, neutrallosses::Bool)
    if neutrallosses
        frags = Iterators.flatten(allneutrallosses.(frags))
    end
    (ionize(frag, c) for frag in frags for c in c_max:-1:1)
end
export alltransitions
function alltransitions(pi::PeptideIon, nterm::ProtNTerm, cterm::ProtCTerm; min_length::Int = 1, precursor::Bool = false, internal::Bool = false, neutrallosses::Bool = false)
    frags = allfragments(pepfragment(pi), nterm, cterm; min_length = min_length, precursor = precursor, internal = internal)
    allions_helper(frags, charge(pi), neutrallosses)
end
function alltransitions(pi::PeptideIon, nterm::ProtNTerm; min_length::Int = 1, precursor::Bool = false, neutrallosses::Bool = false)
    frags = allfragments(pepfragment(pi), nterm; min_length = min_length, precursor = precursor)
    allions_helper(frags, charge(pi), neutrallosses)
end
function alltransitions(pi::PeptideIon, cterm::ProtCTerm; min_length::Int = 1, precursor::Bool = false, neutrallosses::Bool = false)
    frags = allfragments(pepfragment(pi), cterm; min_length = min_length, precursor = precursor)
    allions_helper(frags, charge(pi), neutrallosses)
end

export mass
mass(x::Missing) = 0
function mass(ppos::PeptidePos) mass(ppos.aa) + mass(ppos.label) + mass(ppos.mod) end
function mass(pnterm::PeptideNTerm) mass(pnterm.nterm) + mass(pnterm.mod) end
function mass(pcterm::PeptideCTerm) mass(pcterm.cterm) + mass(pcterm.mod) end
function mass(p::PeptideBase) +(mass.(sequence(p))..., mass(nterm(p)), mass(cterm(p))) end
function mass(p::PeptideFragment) +(invoke(mass, Tuple{PeptideBase}, p), mass.(neutrallosses(p))...) end
function mass(pi::PeptideIon) mass(pepfragment(pi)) + mass(pi.adduct) * pi.adduct_count end

export avgmass
avgmass(x::Missing) = 0
function avgmass(ppos::PeptidePos) avgmass(ppos.aa) + avgmass(ppos.label) + avgmass(ppos.mod) end
function avgmass(pnterm::PeptideNTerm) avgmass(pnterm.nterm) + avgmass(pnterm.mod) end
function avgmass(pcterm::PeptideCTerm) avgmass(pcterm.cterm) + avgmass(pcterm.mod) end
function avgmass(p::PeptideBase) +(avgmass.(sequence(p))..., avgmass(nterm(p)), avgmass(cterm(p))) end
function avgmass(p::PeptideFragment) +(invoke(avgmass, Tuple{PeptideBase}, p), avgmass.(neutrallosses(p))...) end
function avgmass(pi::PeptideIon) avgmass(pepfragment(pi)) + avgmass(pi.adduct) * pi.adduct_count end

export formula
formula(::Missing) = zero(Dict{ProtElement,Int64})
# function formula(pi::PeptideIon) +(formula.(sequence(pi))..., formula(pi.adduct) * pi.adduct_count) end
# function formula(p::PeptideBase) +(formula.(sequence(p))..., formula(nterm(p)), formula(cterm(p))) end
# function formula(ppos::PeptidePos) +(formula(ppos.aa), formula.(ppos.mods)...) end
# function formula(nl::UnimodNeutralLoss) nl.formula end
# function formula(pnterm::PeptideNTerm) ismissing(pnterm.mod) ? formula(pnterm.nterm) : formula(pnterm.nterm) + formula(pnterm.mod) end
# function formula(pcterm::PeptideCTerm) ismissing(pcterm.mod) ? formula(pcterm.cterm) : formula(pcterm.cterm) + formula(pcterm.mod) end
function formula(ppos::PeptidePos) sum(formulas(ppos)) end
function formula(pnterm::PeptideNTerm) sum(formulas(pcterm)) end
function formula(pcterm::PeptideCTerm) sum(formulas(pcterm)) end
function formula(p::PeptideBase) sum(formulas(p)) end
function formula(pi::PeptideIon) sum(formulas(pi)) end

function formulas(ppos::PeptidePos) vcat(formula(ppos.aa), formula(ppos.label), formula(ppos.mod)) end
function formulas(pnterm::PeptideNTerm) vcat(formula(pnterm.nterm), formula(pnterm.mod)) end
function formulas(pcterm::PeptideCTerm) vcat(formula(pcterm.cterm), formula(pcterm.mod)) end
function formulas(p::PeptideBase) vcat(formulas.(sequence(p))..., formulas(nterm(p)), formulas(cterm(p))) end
function formulas(p::PeptideFragment) vcat(invoke(formulas, Tuple{PeptideBase}, p), formula.(neutrallosses(p))) end
function formulas(p::PeptideFragmentNTerm) vcat(formulas.(sequence(p))..., formulas(nterm(p)), formula(cterm(p)), formula.(neutrallosses(p))) end
function formulas(p::PeptideFragmentCTerm) vcat(formulas.(sequence(p))..., formula(nterm(p)), formulas(cterm(p)), formula.(neutrallosses(p))) end
function formulas(p::PeptideFragmentInternal) vcat(formulas.(sequence(p))..., formula(nterm(p)), formula(cterm(p)), formula.(neutrallosses(p))) end
function formulas(pi::PeptideIon) vcat(formulas(pepfragment(pi)), formula(pi.adduct) * pi.adduct_count) end

function parse(::Type{PrecursorIon}, str::AbstractString)
    str_split = match(r"^(.*?)([+-]+)$", str)
    peptide = parse(Peptide, str_split.captures[1])
    charge = str_to_charge(str_split.captures[2])
    ionize(peptide, charge)
end
function parse(::Type{Peptide}, str::AbstractString)
    match_n_seq_c = match(r"^(\w+|(?:\[.*?[^\\]\]))-(.*)-(\w+|(?:\[.*?[^\\]\]))$", str)
    if !isnothing(match_n_seq_c)
        nterm = match_n_seq_c.captures[1]
        str = match_n_seq_c.captures[2]
        cterm = match_n_seq_c.captures[3]
    else
        nterm = "H"
        cterm = "OH"
    end
    if !occursin(r"^([A-Z](?:\[.*?[^\\]\])?(?:\[.*?[^\\]\])?)+$", str)
        error("Peptide has invalid sequence: ", str)
    end
    pos_strs = getfield.(eachmatch(r"[A-Z](?:\[.*?[^\\]\])?(?:\[.*?[^\\]\])?", str), :match)
    poss = vcat(
        parse(PeptidePos, pos_strs[1], ProtNTerm),
        parse.(PeptidePos, pos_strs[2:(end-1)]),
        parse(PeptidePos, pos_strs[end], ProtCTerm)
    )
    pnterm = parse(PeptideNTerm, nterm)
    pcterm = parse(PeptideCTerm, cterm)

    Peptide(poss, pnterm, pcterm)
end
function ppos_parser(str::AbstractString)
    aa = parse(ProtAA, str[[1]])
    label::Union{UnimodLabel, Missing} = missing
    mod::Union{Unimod, Missing} = missing

    if length(str) > 1
        (str1, str2) = match(r"(?:\[(.*?[^\\])\])?(?:\[(.*?[^\\])\])?", str[2:end]).captures

        if !isnothing(str2)
            label = parse(UnimodLabel, str1)
            mod = parse(Unimod, str2)
        elseif !isnothing(str1)
            try label = parse(UnimodLabel, str1)
            catch e mod = parse(Unimod, str1)
            end
        else
            error()
        end
    end

    (aa = aa, label = label, mod = mod)
end
function parse(::Type{PeptidePos}, str::AbstractString)
    x = ppos_parser(str)
    specmod = ismissing(x.mod) ? missing : specify(x.mod, x.aa)
    PeptidePos(x.aa, x.label, specmod)
end
function parse(::Type{PeptidePos}, str::AbstractString, t::Type{T}) where T
    x = ppos_parser(str)
    specmod = ismissing(x.mod) ? missing : specify(x.mod, x.aa, t)
    PeptidePos(x.aa, x.label, specmod)
end
function parse(::Type{PeptideNTerm}, str::AbstractString)
    if str[1] == '[' && str[end] == ']'
        mod = parse(Unimod, replace(str[2:(end-1)], "\\]" => "]"))
        specmod = specify(mod, P_AA_N_term, ProtNTerm)
        PeptideNTerm(P_NT_H, specmod)
    else
        PeptideNTerm(parse(ProtNTerm, str[1:(end)]), missing)
    end
end
function parse(::Type{PeptideCTerm}, str::AbstractString)
    if str[1] == '[' && str[end] == ']'
        mod = parse(Unimod, replace(str[2:(end-1)], "\\]" => "]"))
        specmod = specify(mod, P_AA_C_term, ProtCTerm)
        PeptideCTerm(P_CT_OH, specmod)
    else
        PeptideCTerm(parse(ProtCTerm, str), missing)
    end
end

charge_to_str(x::Int) = repeat(x >= 0 ? '+' : '-', abs(x))
str_to_charge(x::AbstractString) = count("+", x) - count("-", x)

export ordinate
function ordinate(p::PeptideBase) length(sequence(p)) end

nldesc(x::AANeutralLoss) = string(x)
nldesc(x::UnimodNeutralLoss) = string(x)
# nldesc(x::AANeutralLoss) = string(round(Int, mass(x)))
# nldesc(x::UnimodNeutralLoss) = string(round(Int, mass(x)))
export fragmenttype
fragmenttype(p::PeptideFragmentFull) = "M"
fragmenttype(p::PeptideFragmentNTerm) = string(cterm(p))
fragmenttype(p::PeptideFragmentCTerm) = string(nterm(p))
fragmenttype(p::PeptideFragmentInternal) = string(nterm(p)) * string(cterm(p))

# TODO check fixes
export fragmentdesc
function fragmentdesc(p::PeptideFragmentFull)
    nls = collect(skipmissing(vcat(nlsequence(p)..., nlnterm(p), nlcterm(p))))
    *(
        fragmenttype(p),
        (isempty(nls) ? "" : sort(nldesc.(nls)))...
    )
end
function fragmentdesc(p::PeptideFragmentNTerm)
    nls = collect(skipmissing(vcat(nlsequence(p)..., nlcterm(p))))
    *(
        fragmenttype(p),
        string(ordinate(p)),
        (isempty(nls) ? "" : sort(nldesc.(nls)))...
    )
end
function fragmentdesc(p::PeptideFragmentCTerm)
    nls = collect(skipmissing(vcat(nlsequence(p)..., nlnterm(p))))
    *(
        fragmenttype(p),
        string(ordinate(p)),
        (isempty(nls) ? "" : sort(nldesc.(nls)))...
    )
end
function fragmentdesc(p::PeptideFragmentInternal)
    nls = collect(skipmissing(nlsequence(p)))
    *(
        fragmenttype(p),
        string(p.fragment_pos_n), ":", string(p.fragment_pos_c - 1),
        (isempty(nls) ? "" : sort(nldesc.(nls)))...
    )
end

export iondesc
iondesc(pi::FragmentIon) = fragmentdesc(pepfragment(pi)) * charge_to_str(charge(pi))

export textrepresentation
function textrepresentation(pi::PrecursorIon; modsasmasses = false)
    sprint(io -> print(io, pi; modsasmasses = modsasmasses))
end
function show(io::IO, pi::PrecursorIon; modsasmasses = false)
    print(io, peptide(pi); modsasmasses = modsasmasses)
    print(io, charge_to_str(charge(pi)))
end
function textrepresentation(pi::PeptideIon; modsasmasses = false)
    sprint(io -> print(io, pi; modsasmasses = modsasmasses))
end
function print(io::IO, pi::PeptideIon; modsasmasses = false)
    print(io, pepfragment(pi); modsasmasses = modsasmasses)
    print(io, charge_to_str(charge(pi)))
end
function textrepresentation(ps::PeptideFragment; modsasmasses = false)
    sprint(io -> print(io, ps; modsasmasses = modsasmasses))
end
function print(io::IO, ps::PeptideFragment; modsasmasses = false)
    p = peptide(ps)
    print(io, p; modsasmasses = modsasmasses)
    print(io, "; ", fragmentdesc(ps))
end
function textrepresentation(p::Peptide; modsasmasses = false)
    sprint(io -> print(io, p; modsasmasses = modsasmasses))
end
function print(io::IO, p::Peptide; modsasmasses = false)
    showterms = nterm(p).nterm != P_NT_H || !ismissing(nterm(p).mod) || cterm(p).cterm != P_CT_OH || !ismissing(cterm(p).mod)
    showterms && print(io, nterm(p); modsasmasses = modsasmasses)
    print(io, sequence(p); modsasmasses = modsasmasses)
    showterms && print(io, cterm(p); modsasmasses = modsasmasses)
end
function textrepresentation(pseq::PeptideSequence; modsasmasses = false)
    sprint(io -> print(io, pseq; modsasmasses = modsasmasses))
end
function print(io::IO, pseq::PeptideSequence; modsasmasses = false)
    for ppos in pseq
        print(io, ppos; modsasmasses = modsasmasses)
    end
end
function textrepresentation(ppos::PeptidePos; modsasmasses = false)
    sprint(io -> print(io, ppos; modsasmasses = modsasmasses))
end
# function show(io::IO, ppos::PeptidePos)
#     show(io, ppos.aa)
#     !ismissing(ppos.label) && print(io, "\n\t[", string(ppos.label), "]")
#     !ismissing(ppos.mod) && print(io, "\n\t[", string(ppos.mod), "]")
# end
function print(io::IO, ppos::PeptidePos; modsasmasses = false)
    print(io, ppos.aa)
    if !ismissing(ppos.label)
        print(io, "[")
        print(io, ppos.label; asmass = modsasmasses)
        print(io, "]")
    end
    if !ismissing(ppos.mod)
        print(io, "[")
        print(io, ppos.mod; asmass = modsasmasses)
        print(io, "]")
    end
end
function textrepresentation(pnterm::PeptideNTerm; modsasmasses = false)
    sprint(io -> print(io, pnterm; modsasmasses = modsasmasses))
end
# function show(io::IO, pnterm::PeptideNTerm)
#     show(io, pnterm.nterm)
#     ismissing(pnterm.mod) || show(io, pnterm.mod)
# end
function print(io::IO, pnterm::PeptideNTerm; modsasmasses = false)
    if ismissing(pnterm.mod)
        print(io, pnterm.nterm)
    else
        pnterm.nterm == P_NT_H || print(io, pnterm.nterm)
        print(io, "[")
        print(io, pnterm.mod; asmass = modsasmasses)
        print(io, "]")
    end
    print(io, "-")
end
function textrepresentation(pcterm::PeptideCTerm; modsasmasses = false)
    sprint(io -> print(io, pcterm; modsasmasses = modsasmasses))
end
# function show(io::IO, pcterm::PeptideCTerm)
#     show(io, pcterm.cterm)
#     ismissing(pcterm.mod) || show(io, pcterm.mod)
# end
function print(io::IO, pcterm::PeptideCTerm; modsasmasses = false)
    print(io, "-")
    if ismissing(pcterm.mod)
        print(io, pcterm.cterm)
    else
        pcterm.cterm == P_CT_OH || print(io, pcterm.cterm)
        print(io, "[")
        print(io, pcterm.mod; asmass = modsasmasses)
        print(io, "]")
    end
end

export setsil
setsil(ppos::PeptidePos, label::UnimodLabel = fullsil(ppos.aa)) = @set ppos.label = label
setsil(p::Peptide, pos::Int) = @set p.sequence[pos] = setsil(p.sequence[pos])
setsil(p::Peptide, pos::Int, label::UnimodLabel) = @set p.sequence[pos] = setsil(p.sequence[pos], label)
function setsil(p::Peptide, pos::AbstractVector{Int})
    seq = copy(p.sequence)
    for i in pos
        seq[i] = setsil(seq[i])
    end
    @set p.sequence = seq
end
setsil(p::Peptide) = setsil(p, keys(sequence(p)))

export clearsil
clearsil(ppos::PeptidePos) = @set ppos.label = missing
clearsil(p::Peptide, pos::Int) = @set p.sequence[pos].label = missing
function clearsil(p::Peptide, pos::AbstractVector{Int})
    seq = copy(p.sequence)
    for i in pos
        seq[i] = clearsil(seq[i])
    end
    @set p.sequence = seq
end
clearsil(p::Peptide) = clearsil(p, keys(p.sequence))
clearsil(p::PeptideFragment, pos::Int) = @set p.peptide.sequence[pos].label = missing
function clearsil(p::PeptideFragment, pos::AbstractVector{Int})
    seq = copy(p.peptide.sequence)
    for i in pos
        seq[i] = clearsil(seq[i])
    end
    @set p.peptide.sequence = seq
end
clearsil(p::PeptideFragment) = clearsil(p, keys(p.peptide.sequence))
clearsil(p::PeptideIon, pos::Int) = @set p.peptide = clearsil(p.peptide, pos)
clearsil(p::PeptideIon, pos::AbstractVector{Int}) = @set p.peptide = clearsil(p.peptide, pos)
clearsil(p::PeptideIon) = @set p.peptide = clearsil(p.peptide)

export setmod
setmod(ppos::PeptidePos, mod::Unimod) = @set ppos.mod = specify(mod, ppos.aa)
setmod(ppos::PeptidePos, mod::Unimod, t::Type{T}) where T = @set ppos.mod = specify(mod, ppos.aa, t)
function setmod(p::Peptide, pos::Int, mod::Unimod)
    pos == 1 &&
        return @set p.sequence[pos] = setmod(p.sequence[pos], mod, ProtNTerm)
    pos == length(p.sequence) &&
        return @set p.sequence[pos] = setmod(p.sequence[pos], mod, ProtCTerm)
    @set p.sequence[pos] = setmod(p.sequence[pos], mod)
end

export clearmod
clearmod(ppos::PeptidePos) = @set ppos.mod = missing
clearmod(p::Peptide, pos::Int) = @set p.sequence[pos].mod = missing
function clearmod(p::Peptide, pos::AbstractVector{Int})
    seq = copy(p.sequence)
    for i in pos
        seq[i] = clearmod(seq[i])
    end
    @set p.sequence = seq
end
# function clearmod(p::Peptide, pos::AbstractVector{Int})
#     seq = copy(p.sequence)
#     for i in pos
#         seq[i] = clearmod(seq[i])
#     end
#     @set p.sequence = seq
# end
clearmod(p::Peptide) = clearmod(p, keys(p.sequence))
clearmod(p::PeptideFragment, pos::Int) = @set p.peptide.sequence[pos].mod = missing
function clearmod(p::PeptideFragment, pos::AbstractVector{Int})
    seq = copy(p.peptide.sequence)
    for i in pos
        seq[i] = clearmod(seq[i])
    end
    @set p.peptide.sequence = seq
end
clearmod(p::PeptideFragment) = clearmod(p, keys(p.peptide.sequence))
clearmod(p::PeptideIon, pos::Int) = @set p.peptide = clearmod(p.peptide, pos)
clearmod(p::PeptideIon, pos::AbstractVector{Int}) = @set p.peptide = clearmod(p.peptide, pos)
clearmod(p::PeptideIon) = @set p.peptide = clearmod(p.peptide)

export bare
bare(p::PeptideBase) = clearsil(clearmod(p))

export applymods
applymods(p::Peptide, mod::SpecifiedUnimod) = applymods(p, [mod])
function applymods(p::Peptide, mods::AbstractVector{SpecifiedUnimod})
    isempty(mods) && return p

    seq = copy(sequence(p))
    pposs_avail = findall(ppos -> ismissing(ppos.mod), seq)

    mods_seq = filter(mod -> mod.specificity.site ∉ (P_AA_N_term, P_AA_C_term), mods)
    mods_nterm = filter(mod -> mod.specificity.site == P_AA_N_term, mods)
    length(mods_nterm) <= 1 || error()
    mods_cterm = filter(mod -> mod.specificity.site == P_AA_C_term, mods)
    length(mods_cterm) <= 1 || error()

    for i in pposs_avail
        for mod in mods_seq
            ppos = seq[i]
            if ppos.aa == mod.specificity.site
                seq[i] = @set ppos.mod = mod
                # only the first matching mod is applied at any one position
                break
            end
        end
    end

    p = @set p.sequence = seq

    for mod in mods_nterm 
        p = @set p.nterm.mod = mod
    end
    for mod in mods_cterm 
        p = @set p.cterm.mod = mod
    end

    p
end

export allvariablemods
allvariablemods(p::Peptide, mod::SpecifiedUnimod) = allvariablemods(p, [mod])
function allvariablemods(p::Peptide, mods::AbstractVector{SpecifiedUnimod})
    seq = sequence(p)
    pposs_avail = findall(ppos -> ismissing(ppos.mod), seq)

    mods_seq = filter(mod -> mod.specificity.site ∉ (P_AA_N_term, P_AA_C_term), mods)
    mods_nterm = filter(mod -> mod.specificity.site == P_AA_N_term, mods)
    mods_cterm = filter(mod -> mod.specificity.site == P_AA_C_term, mods)

    mod_poss = [i => mod for mod in mods for i in pposs_avail if seq[i].aa == mod.specificity.site]
    seqs = Vector{PeptideSequence}()
    for mod_poss in powerset(mod_poss)
        allunique(first.(mod_poss)) || continue
        seq_tmp = copy(seq)
        for (i, mod) in mod_poss
            ppos = seq_tmp[i]
            ppos = @set ppos.mod = mod
            seq_tmp[i] = ppos
        end
        push!(seqs, seq_tmp)
    end
    
    ps = [@set p.sequence = seq for seq in seqs]

    ps_ref = copy(ps)
    for mod in mods_nterm for p in ps_ref
        push!(ps, @set p.nterm.mod = mod)
    end end

    ps_ref = copy(ps)
    for mod in mods_cterm for p in ps_ref
        push!(ps, @set p.cterm.mod = mod)
    end end

    ps
end

@setup_workload begin
    @compile_workload begin
        pep = parse(Peptide, "SM[Oxidation]SAAPPPPPR")
        ion = ionize(pep, 2)
        mass(pep)
        formula(pep)
        mass(ion)
        formula(ion)
        allfragments(pep, P_NT_y, P_CT_b, min_length = 1, precursor = true, internal = true)
        textrepresentation.(alltransitions(ion, P_NT_y, P_CT_b, min_length = 2, precursor = true, neutrallosses = true))
        pn = clearmod(pep, 2)
        pn = setmod(pn, 2, parse(Peptidomics.Unimod, "Oxidation"))
        clearmod(pn)
    end
end

end # module
