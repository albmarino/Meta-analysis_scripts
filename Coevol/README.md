## Traits correlation under phylogenetic non-independence with Coevol

Scripts to select GC3-poor and GC3-rich genes and to run Coevol.

#### Full geneset stats and filtering
`genespersp.sh` calculates the number of markers for each species.<br>
`gcstat.sh` reports GC3 level statistics per gene overall and at clade level.<br>
`gc3plots.R` plots GC3 levels per clade and outputs a list of species to remove.<br>

#### Extract genes post-filtering
`gcstat_filtgenespersp.sh` reports GC3 level statistics on the filtered dataset after species removal.<br>
`genespersp_filtgenespersp.sh` calculates the number of markers for each species on the filtered dataset.<br>
`gc3plots_gcrich.R` outputs lists of the top GC3-rich genes for every clade.<br>
`gc3plots.R` outputs lists of the top 50 GC3-poor genes for every clade.<br>

#### Run Coevol and read output
`concatenate_4coevol.sh` concatenates listed genes.<br>
`make_coevol_table.R` produces a matrix in input format for Coevol.<br>
`run_coevol.sh` and `read_coevol.sh` run Coevol and check output.<br>
`dnds_coevol.R` extracts Coevol dS and dN/dS.<br>


## References

Lartillot, N., & Poujol, R. (2011). A Phylogenetic Model for Investigating Correlated Evolution of Substitution Rates and Continuous Phenotypic Characters. Molecular Biology and Evolution, 28(1), 729–744. doi: 10.1093/molbev/msq244
