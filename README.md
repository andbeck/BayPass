# BayPass Details

This repo contains scripts for applying Gautier's 2015 BayPass algorithm to an 8 population dataset of pooled re-sequenced data.

All Data are now stored on TTYLMF and BigOne

####Master Analysis of Raw Data has the following detail

* Uses ALLELEFILE, SAMPLEFILE, and prENVFILE (principal components ENVFILE).
* Has outprefix of `ana_prENVFILE`
* run from g_baypass via `BayPass_BaseModel.txt` file

####Simulated Data are Made from this

* uses `simulate.baypass()` function provided by Gautier
* produced a `GPool.pod` file to re-run using g_baypass (must run again to save 3.1.17)
* output is `anaPOD_`  files

####Analysis proceeds by comparing real and simulated data output

* found in `baypass_daph_interp.R` by Becks

####Jacknifed Data created from master, dropping a pond

* used `baypass_jack_simDat.R` file by Becks
* Created ALLELFILE and prENVFILE with `_no<PONDNAME>` suffixes


####Re-analysed all Jacknifed Data

* using again g_baypass via `BayPass_BaseModel.txt` file
* output are `ana.no<PONDNAME>_` files
* did this twice (v2)

#### Simulated from all jacknifed data

* using 'baypass_jack_simDat.R` again
* creates `Gpool.no_<POND>_pods` files
* all other gubbins stored in `JackSimDatOuts` folder

### Major RData files

* `datUse.clean.Rdata` = raw data going in on alleles
* `BayPassInterpret.RData` = outputs from initial analysis of raw data
* `outliers.RData` = outliers tables from inital anlysis of raw data
* `JackSimDat.RData` = all simulated jack datas

### NOTE on parallel

quite easy to implement g_baypass in terminal using parallel via

`parallel < ~/GitHubRepos/BayPass/sevenJack.txt`

and

`parallel < ~/GitHubRepos/BayPass/sevenJackSim.txt`

where script (e.g. sevenJack.txt) is a simple text file of the n (7) calls you want parallelised and as long as paths are specified, or you are in the path where you want it all to happen and where all the files are.....

