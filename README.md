# CandidateVectorSearch_template

This repository serves as a template for implementing a custom computational
backend for [*CandidateSearch*](https://github.com/hgb-bin-proteomics/CandidateSearch)
and the non-cleavable crosslink search in [*MS Annika*](https://github.com/hgb-bin-proteomics/MSAnnika).

An exemplary custom implementation is given in *findTopCandidates* (in `VectorSearch/dllmain.cpp`).
A demonstration of this example implementation is given in `spmv_example.md`.

Users do not need to implement every function in `VectorSearch/dllmain.cpp` but
simply choose the one which function signature fits best to their use case. The
implemented function can then be called from *CandidateSearch* or *MS Annika*
using the parameter value that maps to that function (given in the function's
documentation).

It might be useful to also look at the official implementation of
[*CandidateVectorSearch*](https://github.com/hgb-bin-proteomics/CandidateVectorSearch)
which uses [Eigen](https://eigen.tuxfamily.org/)
for computing matrix products.

## Potential new backends

A (non-exhaustive) list of potential new backends that could be implemented is
given in [Issues](https://github.com/hgb-bin-proteomics/CandidateVectorSearch_template/issues).

## Citing

If you are using [parts of] *CandidateVectorSearch* or this template please cite this [publication](https://doi.org/10.1038/s42004-024-01386-x):

```
Proteome-wide non-cleavable crosslink identification with MS Annika 3.0 reveals the structure of the C. elegans Box C/D complex
Micha J. Birklbauer, Fränze Müller, Sowmya Sivakumar Geetha, Manuel Matzinger, Karl Mechtler, and Viktoria Dorfer
Communications Chemistry 2024 7 (300)
DOI: 10.1038/s42004-024-01386-x
```

## License

- [MIT](https://github.com/hgb-bin-proteomics/CandidateVectorSearch_template/blob/master/LICENSE)

## Contact

[micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
