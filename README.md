# `recom-sim`: A hybrid simulation tool

## Uses

`recom-sim` is a tool based on Nielsen et al.'s `HYBRIDLAB` for simulates hybrids from the genetic data of two reference populations using allele frequencies. It allows for faster simulations, significantly more markers, and the creation of common introgression classes (B1, B2, F2) in one step. While mainly geared towards SNP data, any biallelic data in `GENEPOP` format can be used.

### Example use

```python recom-sim.py input_file.txt 1 --num-offs 1000 --out output_file.txt```

## Input file

Similar to `HYBRIDLAB`, `recom-sim` uses the `GENEPOP` format as input. Any variation of comma or newline separated loci, three or two digit alleles, and tab or space separated data is accepted. Ensure that only two populations are present in the file.

## Options
```
Usage: recom-sim [OPTIONS]

Positional
  input_file          file path to GENEPOP reference
  introgression       choice of {1,2,3}

Optional
  --num_offs INT      number of offspring to simulation (def = 100)
  --p1name   TEXT     name for pop1, used to differentiate backcrosses (def = POP1)
  --p2name   TEXT     name for pop2, used to differentiate backcrosses (def = POP2)
  --exclude           excludes reference populations from
  --out      TEXT     output file name (def = 'out')
  ```

### Introgression level

The positional argument 'introgression' signifies the generational level to simulate.

1. Only F1 hybrids
2. F1, F2, backcrosses to both reference populations: B1POP1, B1POP2 (F1 x POP1, F1 x POP2)
3. F1, F2, B1POP1, B1POP2, B2POP1 and B2POP2 (B1POP1 x POP1, B1POP2 x POP2)

Other introgression classes can be simulated by editing the output and feeding it back into the program (e.g. inputing two 'populations' of F2 with an introgression level of 1 to get F3)

## References
Nielsen, E. E., L. A. Bach and P. Kotlicki (2006). "HYBRIDLAB (version 1.0): a program for generating simulated hybrids from population samples." Molecular Ecology Notes 6(4): 971-973
