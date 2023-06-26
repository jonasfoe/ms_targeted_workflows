module MSUtils2

using DataFrames, Peptidomics, PrecompileTools
import CSV
import Combinatorics.powerset

include("utils.jl")
include("peptides.jl")
include("external_tools.jl")
include("peptide_data.jl")

@setup_workload begin
    pep = parse(Peptide, "SM[Oxidation]SAAPPPPPR")
    ion = ionize(pep, 2)
    @compile_workload begin
        seq_from_skyline_to_proteomics(seq_from_proteomics_to_skyline("SM[Oxidation]SAAPPPPPR"))
        transition_df_from_precursors([ion])
        merge_dicts_to_df(a = Dict(:a => 1), b = Dict(:b => 1))
    end
end

# function is_skyline_fragment(ions::AbstractVector{Ion}, include_precursor::Bool = false)
#     ion_types = "abcxyz"
#
#     if include_precursor
#         ion_types = "M" * ion_types
#     end
#
#     for ion in ions
#         if (ion.ion_type in "Mabcxyz") & (ion.ion_mod in ["", "-H20", "-NH3"])
#             return true
#         end
#     end
#
#     return false
# end
# function is_skyline_fragment(annotation::AbstractString)
#     is_skyline_fragment(split(annotation, " | "))
# end
# function is_skyline_fragment(fragments::AbstractVector{<: AbstractString}, include_precursor::Bool = false)
#     if include_precursor
#         check_true = [
#             r"^[Mabcxyz]\d+\++$",
#             r"^[Mabcxyz]\d+-H2O\w*\++$",
#             r"^[Mabcxyz]\d+-NH3\w*\++$",
#         ]
#     else
#         check_true = [
#             r"^[abcxyz]\d+\++$",
#             r"^[abcxyz]\d+-H2O\w*\++$",
#             r"^[abcxyz]\d+-NH3\w*\++$",
#         ]
#     end
#
#     for fragment in fragments
#         for s in check_true
#             if occursin(s, fragment)
#                 return true
#             end
#         end
#     end
#
#     return false
# end

end # module
