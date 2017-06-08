# Codons
# ======
#
# Compact codon sequence.
#
# A Codon is a short 3 nucleotide sequence, without any 'N's, packed in a
# single 8 bit value.
#
# While BioSequence is an efficient general-purpose sequence
# representation, Codon is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

bitstype 8 Codon{T<:NucleicAcid} <: Sequence
typealias DNACodon Codon{DNA}
typealias RNACodon Codon{RNA}

# Conversion from/to integers
# ---------------------------

Base.convert{T<:NucleicAcid}(::Type{Codon{T}}, nt::UInt8) = reinterpret(Codon{T}, nt)
Base.convert(::Type{UInt8}, cdn::Codon) = reinterpret(UInt8, cdn)
Base.convert{T<:Number}(::Type{T}, cdn::Codon) = convert(T, UInt8(cdn))

# Basic Functions
# ---------------

alphabet(::Type{DNACodon}) = (DNA_A, DNA_C, DNA_G, DNA_T)
alphabet(::Type{RNACodon}) = (RNA_A, RNA_C, RNA_G, RNA_U)
Base.length(cdn::Codon) = 3
Base.eltype{T}(::Type{Codon{T}}) = T
@inline function inbounds_getindex{T}(x::Codon{T}, i::Integer)
    return reinterpret(T, 0x01 << ((UInt8(x) >> (2 + 2(3 - i))) & 0b11))
end
Base.summary{T<:NucleicAcid}(x::Codon{T}) = string(T, " Codon")
Base.:-{T<:NucleicAcid}(x::Codon{T}, y::Integer) = Codon{T}(UInt8(x) - y % UInt8)
Base.:+{T<:NucleicAcid}(x::Codon{T}, y::Integer) = Codon{T}(UInt8(x) + y % UInt8)
Base.:+{T<:NucleicAcid}(x::Integer, y::Codon{T}) = y + x
Base.:(==){T<:NucleicAcid}(x::Codon{T}, y::Codon{T}) = UInt8(x) == UInt8(y)
Base.isless{T<:NucleicAcid}(x::Codon{T}, y::Codon{T}) = isless(UInt8(x), UInt8(y))
