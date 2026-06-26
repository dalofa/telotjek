# Telotjek
Should you inspect the ends of your Streptomycete assembly? If the linear replicons do not contain any of the telomere clusters described in Faurdal et al. 2025 [1] or the six classes of telomeres described in Algora-Gallardo et al 2021 [2] it might be a good idea.

Telotjek provides a quick and dirty way to approach this check. It runs a BLAST against a representative of the 129 telomere clusters generated in Faurdal et al. 2025 [1] or a representative of the six classes outlined in Algora-Gallardo et al 2021 [2] and reports if it matches any of the classes and if it is present in full or truncated. 

If you do not have a match in your *Complete* assembly consider checking if the reads support it being extended further [1] or if you have a potentially novel telomere class on your hands.

## Usage
```
telotjek.py assembly.fasta
```

To use the 6 representatives outlined by Algora-Gallardo et al 2021:
```
telotjek.py assembly.fasta --telo_class
```

## Dependencies
Requires blast

## Source:
1. Faurdal D, Booth TJ, Weber T, Jørgensen TS. Tying up loose ends: recovering thousands of missing telomeres from Streptomyces and other Streptomycetaceae genomes. bioRxiv [Preprint]. 2025 Oct 14:2025.10.14.682034. doi:10.1101/2025.10.14.682034.
2. Algora-Gallardo L, Schniete JK, Mark DR, Hunter IS, Herron PR. Bilateral symmetry of linear streptomycete chromosomes. Microb Genom. 2021 Nov;7(11):000692. doi: 10.1099/mgen.0.000692. PMID: 34779763; PMCID: PMC8743542.
