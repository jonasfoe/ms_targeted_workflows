function is_irt_peptide(pep::Peptide)
	clearsil(pep) in getindex.(irt_peptides, 2)
end

# irts::AbstractVector{<: Tuple{Real, AbstractString, Vararg{AbstractString}}}
function gen_irt_base(irts = irt_peptides)
	# To make sure the iRT values are at least similar to the given scale by the vendor,
	# we take the first and last given irts and anchor their values to the vendor ones.
    irts = sort(collect(irts))

	# find the indices of the first and last given irts in the vendor list
    start_irt = findfirst(x -> x[2] == irts[1][2], irt_peptides)
    end_irt = findfirst(x -> x[2] == irts[end][2], irt_peptides)

    irts_times = first.(irts)
    irts_times .= irts_times .- irts_times[1]
    irts_times .= irts_times ./ irts_times[end]
    irts_times .= irts_times .* (irt_peptides[end_irt][1] - irt_peptides[start_irt][1]) .+ irt_peptides[start_irt][1]

    return zip(irts_times, getindex.(irts, 2))
end

# irts::AbstractVector{<: Tuple{Real, AbstractString, Vararg{AbstractString}}}
# irt_base::AbstractVector{<: Tuple{Real, AbstractString, Vararg{AbstractString}}}
function gen_irt_transformer(irts = irt_peptides, irt_base = gen_irt_base(irt_peptides))
    irts = sort(collect(irts))
    irt_base = sort(collect(irt_base))

    ids = indexin(getindex.(irts, 2), getindex.(irt_base, 2))

	#  for both we make sure to exclude missing irts
	irts = irts[ids .!== nothing]
	irt_base = irt_base[filter(id -> id !== nothing, ids)]

	@assert length(irts) == length(irt_base) >= 2

    pre::Vector{Float64} = first.(irts)
    post::Vector{Float64} = first.(irt_base)
    return let pre = pre, post = post
        function(x::Float64)
            i = findfirst(pre_val -> x < pre_val, pre)

            if i == nothing
                i = lastindex(pre)
            end

			if i > 1
            	i -= 1
			end

            return post[i] + ((x - pre[i]) / (pre[i + 1] - pre[i])) * (post[i + 1] - post[i])
        end
    end
end

skyline_translations = [
	# "M[+240.147393]" => "tmt-oxM",
	# "K[+232.166677]" => "tmt-heavyK",
	# "R[+234.160747]" => "tmt-heavyR",
	# "I[+231.169642]" => "tmt-heavyI",
	# "L[+231.169642]" => "tmt-heavyL",
    # "P[+230.166287]" => "tmt-heavyP",
    # "V[+230.166287]" => "tmt-heavyV",
    # "K[+232.166677]" => "tmtheavyK",
    "M[+15.994915]"  => "M[Oxidation]",
    # "K[+43.005814]"  => "carbaK",
    "C[+57.021464]"  => "C[Carbamidomethyl]",
    # "C[+324.216141]" => "itmtC",
    # "K[+224.152478]" => "tmtK",
    # "H[+224.152478]" => "tmtH",
    # "S[+224.152478]" => "tmtS",
    # "T[+224.152478]" => "tmtT",
    "K[+8.014199]"   => "K[Label:13C(6)15N(2)]",
    "R[+10.008269]"  => "R[Label:13C(6)15N(4)]",
    "I[+7.017164]"   => "I[Label:13C(6)15N(1)]",
    "L[+7.017164]"   => "L[Label:13C(6)15N(1)]",
    "F[+10.027228]"   => "F[Label:13C(9)15N(1)]",
    "P[+6.013809]"   => "P[Label:13C(5)15N(1)]",
    "V[+6.013809]"   => "V[Label:13C(5)15N(1)]",
    "A[+4.007099]"   => "A[Label:13C(3)15N(1)]",
    # "[+224.152478]"  => "tmt-",
    # "[+43.005814]"  => "carba-","A[]" = "heavyA"
]

function seq_from_skyline_to_proteomics(sequence::AbstractString)
    for pair in skyline_translations
		# mods with unspecific AA
		if startswith(pair[1], '[') & endswith(pair[1], ']')
			# c-terminal
			if startswith(pair[2], '-')
				if endswith(sequence, pair[1])
					sequence = sequence[1:(end - length(pair[1]))] * pair[2]
				end
			else
				seq_v = split(sequence, pair[1])
				if length(seq_v) > 1
					sequence = join(x == "" ? "" : x[1:(end - 1)] * pair[2] * x[end] for x in seq_v)
				end
			end
		else
			sequence = replace(sequence, pair)
		end
    end

    return sequence
end

function seq_from_proteomics_to_skyline(sequence::AbstractString)
    for pair in skyline_translations
		# mods with unspecific AA
		if startswith(pair[1], '[') & endswith(pair[1], ']')
			# c-terminal
			if startswith(pair[2], '-')
				if endswith(sequence, pair[2])
					sequence = sequence[1:(end - length(pair[2]))] * pair[1]
				end
			else
				seq_v = split(sequence, pair[2])
				if length(seq_v) > 1
					sequence = join(x == "" ? "" : x[1] * pair[1] * x[2:end] for x in seq_v)
				end
			end
		else
			sequence = replace(sequence, reverse(pair))
		end
	end

    return sequence
end
