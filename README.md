# IDEPA-XMBD

<div align=center><img src="./figs/IDEPA_figs.png" width="50%" height="50%" ></div>
&nbsp;


We evaluated five state-of-the-art tools (RankComp v1/v2, PenDA, Peng, and Quantile) through classic computational (precision, Type one error control, parameter evaluation, robustness, and similarity) and functional (pathway enrichment and survival analysis) criteria. We also integrated these tools into a user-friendly tool kit, IDEPA-XMBD , to facilitate individualized DEAs in proteomics studies.

A pre-print describing the method is available at bioRxiv: [Application of personalized differential expression analysis in human cancer proteome](https://www.biorxiv.org/content/10.1101/2021.07.18.452812v2)

## Install
We use docker to encapsulate the command line version and the plotly version of IDEPA-XMBD separately.

#### Cmd Version
Pull the docker image of IDEPA-XMBD:
```shell
docker pull ychlouie/idepa_cmd:latest
```

Create a docker container containing IDEPA-XMBD:
```shell
docker run -it ychlouie/idepa_cmd:latest
```

#### Plotly Version
Pull the image:
```shell
docker pull lylan/idepa:latest
```

Create a docker container
```shell

```

## Usage
#### Cmd versioin
After entering the container, use `-h` to view the IDEPA-XMBD module information:
```shell
python /IDEPA-XMBD/individual_depa.py -h

# View specific module information 
python /IDEPA-XMBD/individual_depa.py [module name] -h
```

We also provide sample data for each module ：
```shell
python /IDEPA-XMBD/individual_depa.py [module name] -p /IDEPA-XMBD/parameters_file/test_parameters.txt
```

### Plotly version
After entering the container, you can operate according to [procedure.docx](./procedure.docx)
