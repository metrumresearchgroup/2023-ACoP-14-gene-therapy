# Summary

This is an open science project that aimed to provide open access to a QSP model for hematopoietic stem cell (HSC) differentiation. The model was implemented in both R and julia.

**Goal of this work**: to develop a versatile platform model for HSC differentiation to erythrocytes, lymphocytes, and granulocytes. The model aimed to predict the dynamics of cell reconstitution after ex-vivo gene therapy. 

**Motivation of this work**: there have been many models developed for HSC differentiation into subset of cells, or development of HSC-derived cells (e.g. T cell development in thymus). However, these models failed to capture the complexity and the multi-organ nature of HSC differentiation. Our work aims to fill this gap. The result is an integrated multi-organ QSP model that captured HSC differentiation to red blood cells (RBCs), granulocytes, T cells, and B cells and their dynamics in blood and in secondary lymphoid organs. 

**Method**: the model was built in a stepwise fashion by integrating published models while incorporating multiple layers of progenitors described in literature and calibrate the model based on experimental data obtained in mouse. The model was then scaled to human. 

**Results**: we demonstrate that our model can capture the granulocytes and lymphocytes reconstitution after ex-vivo gene therapy on a patient that suffered from adenosine deaminase‐deficient severe combined immunodeficiency (ADA-SCID), as well as blood cell and hemoglobin reconstitution after a patient with sickle cell disease (SCD) receiving ex-vivo gene therapy. We also demonstrate that we can capture the vector copy number (VCN) dynamics on ADA-SCID mice after receiving ex-vivo gene therapy.

**Future work**: we believe our model can be adapted to predict the outcome under different gene therapy protocol for the following diseases:

- Transfusion-dependent <img src="https://latex.codecogs.com/svg.image?\beta" title=" " />-thalassemia (currently approved in the U.S. under the name [ZYNTEGLO](https://investor.bluebirdbio.com/news-releases/news-release-details/bluebird-bio-announces-fda-approval-zynteglor-first-gene-therapy))
- Pyruvate Kinase Deficiency (PKD) ([NCT04105166](https://clinicaltrials.gov/ct2/show/NCT04105166))
- Metachromatic leukodystrophy (MLD) (currently authorized in Europe under the name [Libmeldy](https://www.ema.europa.eu/en/medicines/human/EPAR/libmeldy))

We also believe our model can benefit from additional calibration/ validation with additional data at progenitor level. Details are discussed in the [future work](#future-work) and [Challenges and uncertainties](#challenges) sections.

**How to use this repo**: see [Content](#content) and [Setup](#r-setup-section) sections.

## Model development

This model was built on 4 published models: 

- A model on human HSC differentiation to red blood cells (RBCs) described in [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638).

- A mouse B cell development model described in [Shahaf et al., 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00077/full).

- A mouse T cell development model described in [Thomas-Vaslin et al., 2008](https://journals.aai.org/jimmunol/article/180/4/2240/84513/Comprehensive-Assessment-and-Mathematical-Modeling).

- A human naive T cell dynamics model described in [Sove et al., 2020](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12546)

<img src="img/Diagram_integration.png" alt="model structure" width="500" >

The model was developed in 4 parts: 

1. Implementation of existing human HSC -> RBC differentiation model described in [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638). <a id='human-scd'></a>

2. Scale the HSC -> RBC differentiation model from human to mouse. <a id='mouse-scd'></a>

3. Integrate a T cells development model from [Thomas-Vaslin et al., 2008](https://journals.aai.org/jimmunol/article/180/4/2240/84513/Comprehensive-Assessment-and-Mathematical-Modeling), a B cells development model from [Shahaf et al., 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00077/full), and HSC -> granulocyte differentiation model into the mouse HSC -> RBC differentiation model in mouse. <a id='mouse-integrated'></a>

4. Scale the integrated mouse HSC differentiation model to human. <a id='human-integrated'></a>

Eventually, our model was able to predict the immune cell reconstitution of a patient with ADA-SCID after receiving ex-vivo gene therapy ([Aiuti et al., 2009](https://www.nejm.org/doi/pdf/10.1056/NEJMoa0805817)). 

In addition, our model was validated using the RBC recovery after ex vivo gene therapy in a patient ([Ribeil et al., 2017](https://www.nejm.org/doi/full/10.1056/nejmoa1609677)), ex vivo gene therapy in mouse ([Carbonaro et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22833548/)), and stem cell/ progenitor transplants in mouse ([Boyer et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30905737/)).

All models in this folder were implemented in R and solved using `mrgsolve`. See [Setup](#r-setup-section) for more information. 

Here we list some key assumptions that goes into the model: 

- Cell trafficking between bone marrow and other organs were instantaneous. 
- The total number of proliferating long-term hematopoietic stem cells (i.e. LT-HSC) was a fixed number.
- Differentiation of cells was irreversible. 
- Proliferation and differentiation dynamics of progenitor cells were comparable between mouse and human. 
- Reticulocytes and red blood cells did not proliferate. 
- Many the rates, e.g. proliferation rate of cells, globin synthesis rate, were constants. 
- Globin synthesis rates were linearly correspond to gene copy number. 
- In sickle cell models, the vector copy number per genome in a transduced CD34+ cell is 0.5. That is to say, for each transduced cell, there is 1 vector.
- The oxygen saturation on adult hemoglobin that carries sickle cell mutation (HbS) is 68%, healthy/ (with or without T87Q mutation) adult hemoglobin (HbA) is 74%, fetal hemoglobin (HbF) is 88%. Hemoglobin A2 (HbA<sub>2</sub>) carries no significant amount of oxygen in healthy human. 
- T cell dynamics were different between adult mice, young mice, and human. 


# Implementation <a id='human-scd-implementation'></a> and validation of red blood cell differentiation in human 

This model was described in [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638). Here we document our implementation of the model and adjustments we made. 


## Model implementation

Compare to the model described in [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638), we have some major changes.

1. Globin synthesis and assembly implementation. This adjustment was to match the unit in both side of the differential equation and match the fluxes in the system. Here we use <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>-globin synthesis in RET as an example. 

    In Zheng et al., 2021, the dynamics was written as 
    
    <img src="https://latex.codecogs.com/svg.image?\frac{d}{dt}&space;\beta_{RET0}&space;=&space;k_{syn\beta}&space;\cdot&space;RET_0&space;-&space;(k_{deg\beta}&space;&plus;&space;k_{RET2RBC})&space;\cdot&space;\beta_{RET0}&space;-&space;k_{on\alpha\beta}&space;\cdot&space;\alpha_{RET0}&space;\cdot&space;\beta_{RET0}&space;/V_{RET0}&space;&plus;&space;k_{off\alpha\beta}&space;\cdot&space;\alpha\beta_{RET0}" title=""  width="600" />

    
    In this equation, the left side and
    <img src="https://latex.codecogs.com/svg.image?k_{syn\beta}&space;\cdot&space;RET_0&space;" title="k_{syn\beta} \cdot RET_0 " /> represented mass change over time, but 
    <img src="https://latex.codecogs.com/svg.image?k_{on\alpha\beta}&space;\cdot&space;\alpha_{RET0}&space;\cdot&space;\beta_{RET0}&space;/V_{RET0}&space;" title="" /> 
    and 
    <img src="https://latex.codecogs.com/svg.image?k_{off\alpha\beta}&space;\cdot&space;\alpha\beta_{RET0}" title="" /> 
    were concentration changes over time. 

    To correct this mismatch, we changed the equation into 
    <img src="https://latex.codecogs.com/svg.image?\frac{d}{dt}&space;\beta_{RET0}&space;=&space;k_{syn\beta}&space;\cdot&space;RET_0&space;-&space;(k_{deg\beta}&space;&plus;&space;k_{RET2RBC})&space;\cdot&space;\beta_{RET0}&space;-&space;k_{on\alpha\beta}&space;\cdot&space;\alpha_{RET0}&space;\cdot&space;\beta_{RET0}&space;&plus;&space;k_{off\alpha\beta}&space;\cdot&space;\alpha\beta_{RET0}&space;\times&space;V_{RET0}" title="\frac{d}{dt} \beta_{RET0} = k_{syn\beta} \cdot RET_0 - (k_{deg\beta} + k_{RET2RBC}) \cdot \beta_{RET0} - k_{on\alpha\beta} \cdot \alpha_{RET0} \cdot \beta_{RET0} + k_{off\alpha\beta} \cdot \alpha\beta_{RET0} \times V_{RET0}" /> 

    Similar changes were carried out in all the globin differential equations. 


2. Long-term haematopoietic stem cell proliferation rate. This rate was not provided in Zheng et al., 2021. The parameter used here was obtained from [Catlin et al., 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3099568/) (0.057 day<sup>-1</sup>).  

3. Alpha globin synthesis rate. This rate was decreased from <img src="https://latex.codecogs.com/svg.image?1.5&space;\times&space;10^{-6}" title="1.5 \times 10^{-6}" /> nmol.cell<sup>-1</sup>day<sup>-1</sup> to <img src="https://latex.codecogs.com/svg.image?6&space;\times&space;10^{-7}" title="6 \times 10^{-7}" /> nmol.cell<sup>-1</sup>day<sup>-1</sup>. 

    We tuned this parameter because the simulation result based on parameters reported in Zheng et al., 2021 did not fit normal human physiology. Simulated red blood cell in a normal person was ~12 trillion, half of what it was supposed to be, yet hemoglobin concentration inside red blood cell was more than twice of the normal value. 

    To tune this parameter, we used 3 metrics reported in [Pittman, 2011](https://pubmed.ncbi.nlm.nih.gov/21634070/): 
        
        i.   RBC count. In healthy human, it should be ~ 4 million per uL.

        ii.  Hemoglobin concentration in blood. In healthy human, it should be 12 ~ 15.5g/dL. 

        iii. Hemoglobin concentration in RBC. This value should be ~ 330g/L. 

    ![alpha globin synthesis rate tuning](img/scd_ksynalpha.png)

    Overall, when the globin synthesis rate equaled <img src="https://latex.codecogs.com/svg.image?6&space;\times&space;10^{-7}" title="6 \times 10^{-7}" /> nmol.cell<sup>-1</sup>day<sup>-1</sup>, most of the observed metrics fell into the normal range. Thus, in our future implementation, we used the tuned value. 

    In addition, this change was overall in line with observed hemoglobin synthesis rate.
    In [Hibbert et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4457298/), the hemoglobin synthesis rate in rat was described as ~ 0.11-0.42mg.g<sup>-1</sup>.day<sup>-1</sup>. Assuming human had the same synthesis rate, a body of 70kg, and 240 billions of reticulocytes (number estimated based on [Koepke and Koepke, 1986](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-2257.1986.tb00093.x) that number of reticulocyte is 1/100 of number of red blood cells), then human's supposed to have a hemoglobin synthesis rate between <img src="https://latex.codecogs.com/svg.image?5.9&space;\times&space;10^{-7}" title="5.9 \times 10^{-7}" /> and <img src="https://latex.codecogs.com/svg.image?2.2&space;\times&space;10^{-6}" title="2.2 \times 10^{-6}" /> nmol.cell<sup>-1</sup>day<sup>-1</sup>. By the assumption that synthesis rates between <img src="https://latex.codecogs.com/svg.image?\alpha, &space;\beta, &space; \gamma, &space; \delta" title=" " /> globins were 1: 0.5: 0.03: 0.04, then the synthesis rate for <img src="https://latex.codecogs.com/svg.image?\alpha" title=" " />-globin should be in the range of <img src="https://latex.codecogs.com/svg.image?3.8&space;\times&space;10^{-7}" title="3.8 \times 10^{-7}" /> to <img src="https://latex.codecogs.com/svg.image?1.4&space;\times&space;10^{-6}" title="1.4 \times 10^{-6}" /> nmol.cell<sup>-1</sup>day<sup>-1</sup>.
    
    Moreover, in SCD patients, it was reported that Hb synthesis rate averages 0.2g.kg<sup>-1</sup>day<sup>-1</sup> ([Hibbert et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4084935/)). This converted to <img src="https://latex.codecogs.com/svg.image?\alpha" title=" " />-globin synthesis rate of <img src="https://latex.codecogs.com/svg.image?6.8&space;\times&space;10^{-7}" title=" " /> nmol.cell<sup>-1</sup>day<sup>-1</sup>, close to the value we identitified. 
    
    ## Model validation

    We validated our implementation in 2 different scenarios: 
    
    1. Steady state

       A summary table of steady state value was as follows. These values were overall in line with what is observed in human ([Hibbert et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4084935/)).

        |                | RBC (#/uL) | Hb (g/dL) | venous oxygen (mL/dL) | Adult hemoglobin (g/dL) |
        | -----------    | -----------|-----------| --------------------- | ---------------- |
        | Healthy person | 4.2M       | 14.1      |  13.2                 | 12.3 (HbA)       |
        | SCD patients   | 1.8M       | 7.7       |  6.6                  | 6.8 (HbS)        |

        We also presented how the system reached these equilibria, starting with 1 HSC. Similar figures were also presented in  [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638) Figure 4. 

        ![healthy_sickle](img/scd_healthy_sickle.png)


    2. Dynamics after a sickle cell patient goes through lentiviral-induced hematopoietic stem cell transplant, as reported in [Ribeil et al., 2017](https://www.nejm.org/doi/full/10.1056/nejmoa1609677).

    In this case, we compared 4 simulated variables with observed data after the transplant: 
    
    - Total RBC count
    - Hemoglobin (Hb) concentration, both HbS and HbA with T87Q mutation (the transgene carried by lentivirus). 
    - reticulocytes (RET) count
    - Hb concentration in RBC

    The simulation results as well as clinical observation were shown as follows. Overall, the steady states of simulation and patients were similar, though the RBC and Hb recovery dynamics is a little off. This might be due to infusion this patient received prior to receiving ex-vivo gene therapy. These infused RBCs/ Hb are present for at least 3 months after gene therapy, judging by HbA level in patient blood. 
    
    Note that the predicted total RET count (cyan line, top right panel) correspond well to the observed data. This indicated the model probably work well upstream of RET. 
    
    Note that the observed HbS and HbA<sup>T87Q</sup> almost converged at the steady state. This would make sense if any of the following cases were true: 1. the VCN per genome was higher than what is assumed; 2. the conditioning strength was greater than assumed;  3. the DNA segment encodes HbA<sup>T87Q</sup> was more accessible than the DNA fragment that encodes HbS, leading to higher transcription and translation rates of HbA<sup>T87Q</sup> than HbS in transduced cells. 

    ![Ribeil et al., 2017](img/scd_Ribeil2017.png)

# Scale the HSC -> RBC differentiation model from human to mouse

The goal of this step was to obtain a HSC -> RBC differentiation model in mouse, so that we could eventually incorporate granulocytes and lymphocyte dynamics based on published model. We'd like to have a model in mouse because there are more literature data related to HSC differentiation in mouse. 

The graphic summary of the model was as follows. Here, we assumed the feedback impact on CFU-E dynamics and Hb saturation level in mouse were comparable to human.

![](img/mouse_erythrocyte_model.png)

This model was based on the following adjustments:

1.  Adjustments on cells' mean residence time

    - Red blood cell (RBC) lifespan. It was reduced from 120 days to 40 days ([Putten and Croon, 1958](https://ashpublications.org/blood/article/13/8/789/33568/The-Life-Span-of-Red-Cells-in-the-Rat-and-the)).

    - Reticulocyte (RET) maturation time. It was decreased to 2 days from 3 days ([Ney, 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3157046/)). 

    - Colony-forming unit erythroid (CFU-E) cell residence time. It was reduced to 2 days from 7 days ([Palis, 2014](https://pubmed.ncbi.nlm.nih.gov/24478716/#:~:text=Definitive%20RBCs%20in%20mammals%20circulate,nucleated%20cells%20before%20ultimately%20enucleating.)). 


2. Adjustments on cell amplification/ division time

    - Short-term hematooietic stem cells (ST-HSC) and multipotent progenitors (MPP) amplification times. They were adjusted to 32 and 280 times ([Busch et al., 2015](https://www.nature.com/articles/nature14242)).

    - Burst-forming unit erythroid (BFU-E) cell amplification time. It was reduced to 16 times ([Palis, 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904103/)).

    - Common myeloid progenitors (CMP) number of division. The amplification was decreased to 8. It was estimated based on the flux from in and out of CMP. 

        - Influx. This is from MPP. Amongst all the efflux from MPP, the ratio between CMP : CLP is 1:180 ([Busch et al., 2015](https://www.nature.com/articles/nature14242)), thus we could assume ~ 114k of MPP goes to CMP per day.  
    
        - Efflux. This was based on the estimation from RBC dynamics. 
        
        Assuming mouse RBC lifespan = 40 days, RBC count = 10M per uL, with 2mL blood, then the total RBC count in mouse = 10M * 2k = 20B. For each day, the RBC death was 20B/ 40 = 500M. Assuming 32 times amplification at CFU-E stage and 16 times of amplification at BFU-E stage, we calculated that CMP efflux to BFU-E ~800k. 

        With influx:efflux ratio ~ 8, we estimated that CMP has 3 times of division (i.e. amplification = 8). 

3. Other modifications

    - All transduced branches were consolidated. This simplification was purely technical and had no scientific impact on the model. 

    - Volumes of RET and RBC.  The scaled down volume of RET/ RBC was 0.05 * 10<sup>-12</sup> L. This number was calculated based on RBC diameter data reported in [Fukuda et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658564/). This number was also in line with what is reported in [Raabe et al., 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3189672/) and [Bollinger and Everds, 2012](https://www.sciencedirect.com/science/article/pii/B9780123820082000143). 


    - <img src="https://latex.codecogs.com/svg.image?\alpha" title="" />-globin synthesis rate. It was reduced to  <img src="https://latex.codecogs.com/svg.image?3.8&space;\times&space;10^{-7}" title="" /> nmol.cell<sup>-1</sup>day<sup>-1</sup>. The assumption was that total Hb synthesis in human, mouse, and rats were comparable; thus, Hb synthesis rate per cell in mouse should be reduced because there were more numerous in mouse. In addition, mouse RBCs were smaller, thus should have less hemoglobin synthesis capacity.

## Model validation

Here we used both steady state value and RBC recovery dynamics after stem cell and progenitor transplant to validate out model. 

Simulated steady state values and their corresponding physiological values in mouse were summarized as follows. Note that all data from [Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4) were adjusted to account for the fact that its BM cells were isolated from tibia and femur, bones that accounted for ~14% of total bone marrow cells ([Colvin et al., 2004](https://www.nature.com/articles/2403268)). We assumed that total bone marrow cell count ~ 500M, and bone marrow cells from tibia and femur was ~60-80M ([Swamydas et  al., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3732092/)). 

| Readout type           |    predicted number  | Ref number |
| -----------     |   --------- |  -----------  |
| RET count  (#/uL) | 479k      | 200 - 500k  ([Bollinger and Everds, 2012](https://www.sciencedirect.com/science/article/pii/B9780123820082000143)) |
| RBC count  (#/ul) |  9.5M      | 10.2M ([Fukuda et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658564/)) |
| Hb conc (g/dL)   |   13.1     | 13.6 - 16.4 ([Raabe et al., 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3189672/)) |
| Hb in RBC  (g/L)  |   274     | 270 - 330  ([Bollinger and Everds, 2012](https://www.sciencedirect.com/science/article/pii/B9780123820082000143)) |
| MPP               |    15k     |  75k-92k ([Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4), [Dong et al., 2019](https://www.mdpi.com/2073-4409/8/9/951/htm))  |
| CMP               |   913k   |  755k- 3M ([Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4), [Dong et al., 2019](https://www.mdpi.com/2073-4409/8/9/951/htm), [Busch et al., 2015](https://www.nature.com/articles/nature14242))  |
| CFU-E amplification time |    32.8 | 16-32 ([Palis, 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904103/)) |   

The discrepancy on the MPP count between simulated and observed values might be due to current model only tracks MPP that will differentiate into the erythroid lineage. These MPPs accounted ~20% of total MPP ([Pietras et al., 2015](https://www.sciencedirect.com/science/article/pii/S193459091500212X), [Wilson et al., 2008](https://www.cell.com/cell/pdf/S0092-8674(08)01386-X.pdf)). 


We then compared the dynamics of cell reconstitution after stem cell and progenitor transplants in mouse. Data was obtained from [Boyer et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30905737/). Briefly, these mice went through sublethal radiation conditioning before being transplanted with a fixed number of GFP-expressing stem cells/ progenitor. 

Here, we assume the sublethal irradiation leads to 10% - 50% loss of progenitor cells in bone marrow. 

![](img/mouse_RBC_HSCT_MPPT.png)

We noticed the steady state of HSC transplant was overall similar to reported data, though the dynamics is slower. This might be due to HSC proliferate faster after a reduction in number [Busch et al., 2015](https://www.nature.com/articles/nature14242). The dynamics after MPP transplant tracked observed data well. 

# Incorporate granulocytes and lymphocytes into the mouse HSC differentiation model 

The main incentive to include lymphoid dynamics is due to our interest in adopt this model for gene therapy for immune deficiency diseases. Granulocytes were included as placeholder so that the MPP count would match observed data. Thus, we took very different approaches between myeloid and lymphoid branches: we built a coarse granulocyte differentiation model follow the differentiation map MPP -> CMP -> GMP -> granulocytes (GM), while developing a more physiologically-based lymphoid branch. 

The graphic summary of the model was provided as follows.

![](img/mouse_full_model.png)


## Incorporate GMP and granulocytes

Here, we only considered 2 layers: granulocyte (GM) and its progenitor cells (GMP) ([Ferrari et al., 2020](https://www.nature.com/articles/s41576-020-00298-5)). GMP was differentiated from CMP. 

Parameter adjustment & selections were listed as follows.

- Differentiation from CMP to GMP. This number was set to be 2.5 <sup>-1</sup>, in the range of 0.04-4 day<sup>-1</sup> reported in [Busch et al., 2015](https://www.nature.com/articles/nature14242). 

- Mean residence time of GMP. This was set to be 0.12 day ([Busch et al., 2015](https://www.nature.com/articles/nature14242)). 

- Turnover of granulocytes. It was indicated that, on average, GM has lifespan of 6-8 hours in mouse, monkey, and human ([Summers et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20620114/), [He et al., 2018](https://www.jimmunol.org/content/200/12/4059), [Galbraith et al., 1965](https://ashpublications.org/blood/article/25/5/683/37897/Patterns-of-Granulocyte-Kinetics-in-Health)). This translated to death rate ~ 3-4 day<sup>-1</sup>. 

- Amplification at ST-HSC and MPP. MPP amplification round was increased from 280 to 450. This number was the upper bound reported in [Busch et al., 2015](https://www.nature.com/articles/nature14242). ST-HSC amplification number was increased to 90. This number was higher than the number reported in [Busch et al., 2015](https://www.nature.com/articles/nature14242). However, this might be due to [Busch et al., 2015](https://www.nature.com/articles/nature14242) and [Zheng et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34139105/) uses different numbers of LT-HSC for parameter estimation.

- Amplification from GMP to GM. This parameter assumed to be 32 because there were 4 additional layers of mitotic progenitors between GMP and neutrophil, the most abundant granulocytes ([Hong, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5662779/)). Assuming each layer of mitotic progenitors came with additional round of cell division, the amplification between GMP & GM would be 2<sup>5</sup> = 32.

- Granulocyte count in mouse blood. This number was reported to ~ 6k per uL ([Nemzek et al., 2001](https://pubmed.ncbi.nlm.nih.gov/11713907/)).Assuming mouse blood volume ~ 2mL, the total number of GM in mouse would be ~ 12 millions. Note that this number may vary hugely between different mouse strings ([von Vietinghoff and Ley, 2008](https://www.jimmunol.org/content/181/8/5183)). 

## Incorporate lymphoid arms

The B cell and T cell development models were taken from literature. To link these literature models with our model, we introduce common lymphoid progenitors (CLP). 


### Incorporate common lymphoid progenitors

To incorporate common lymphoid progenitor (CLP), a progenitor cell that was derived from multipotent progenitor (MPP) cells, we took parameters of differentiation rate, proliferation rate, and death rate from [Busch et al., 2015](https://www.nature.com/articles/nature14242). The parameter of CLP exported to thymus was estimated based on [Zlotoff and Bhandoola, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076003/).


### Incorporate B cell dynamics

The B cell model was taken from [Shahaf et al., 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00077/full). The implementation and documentation of this model is available in [B cell folder](../Bcell/). 

Additional 2 divisions are assumed between CLP and propreB. This number was tuned to match the cell influxes into the propreB compartments that are reported in [Shahaf et al., 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00077/full).

In addition, we used the parameters that is derived from progenitor depleted mice, based on the experimental procedure described in [Busch et al., 2015](https://www.nature.com/articles/nature14242). 

### Incorporate T cell dynamics

The T cell model was taken from [Thomas-Vaslin et al., 2008](https://www.jimmunol.org/content/180/4/2240.long). The implementation and details of verification can be found in [T cell folder](../Tcell/).


The additional 32 times of amplification between CLP and double negative (DN) cells was estimated based on the assumption that there were a round of division when CLP differentiate into early thymic progenitors (ETP), and ETPs divided 8 times in thymus before differentiate into DNs ([Seddon and Yates, 2018](https://pubmed.ncbi.nlm.nih.gov/30129206/)). Note that in our model, we omit ETP for simplicity.


The rate for CLP to be exported to thymus was unknown, presumably because this number was very small thus hard to observe in experimental setting. Thus, we tuned this rate to fit with observed data from 2 different source: 
1. the CLP count that exported to thymus should be ~ 10 cells per day ([Zlotoff and Bhandoola, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076003/));
2. naive T cell count should be less than 1700 cells per uL blood at steady state ([Boyer et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30905737/)); 

A parameter scan suggested that this rate should be <img src="https://latex.codecogs.com/svg.image?2.5\times 10^{-4}" title=" " />day<sup>-1</sup>. 

![](img/mouse_alphaCLP2DNscan.png)


## Model validation

### Steady state

Simulated steady state and their corresponding physiological values were summarized as below. Overall, simulated steady state of the cells was consistent with what is observed in literature. 

| Cell type | Cell count   | Ref number | 
| --------- | ------------ | ---------- |
| MPP       |   60k        |  75k-92k ([Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4), [Dong et al., 2019](https://www.mdpi.com/2073-4409/8/9/951/htm))  |
| CMP       |   2M         | 755k- 3M ([Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4), [Dong et al., 2019](https://www.mdpi.com/2073-4409/8/9/951/htm), [Busch et al., 2015](https://www.nature.com/articles/nature14242))  |
| GMP       |   757k       | ~1M ([Bae et al., 2019](https://www.nature.com/articles/s41467-019-11386-4)), ~3M([Dong et al., 2019](https://www.mdpi.com/2073-4409/8/9/951/htm)) |
| GM        | 8.4k per uL  | 1.95-12.01k per uL ([Nemzek et al., 2001](https://pubmed.ncbi.nlm.nih.gov/11713907/))  |
| RET       | 498k per uL  | 200 - 500k  ([Bollinger and Everds, 2012](https://www.sciencedirect.com/science/article/pii/B9780123820082000143)) |
| RBC       | 9.9M per uL  | 10.2M per uL ([Fukuda et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658564/)) |
| Hb conc   | 13.2g/dL     | 13.6 - 16.4g/dL ([Raabe et al., 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3189672/)) |
| Hb in RBC | 261g/dL      | 270-330g/dL  ([Bollinger and Everds, 2012](https://www.sciencedirect.com/science/article/pii/B9780123820082000143)) |
| lymphocyte | 1.8k     | 0.12 - 24k ([JHU mouse facility](https://researchanimalresources.jhu.edu/wp-content/uploads/2017/09/2017JHPhenoCoreCBCmice.pdf)) |

### Dynamics after HSCT

The model was also validated by looking at the reconsititution dynamics of RBC, B cells, and T cells after HSC transplant in mouse (data obtained from [Busch et al., 2015](https://www.nature.com/articles/nature14242)).

![dynamics_validation](img/mouse_RBC_T_B_HSCT.png)

Overall, the RBC reconstitution was slower, as we observed in our previous simulation. The predicted B cell reconstitution was in the standard deviation of the observed data. 
T cell recovery at a 50% conditioning strength (i.e. loss of 50% progenitors) appeared to match the early reconstitution dynamics the best, but a 30% conditioning yielded closer steady state. There are a couple of possible explanations for the discrepancy: 

1. Difference between mice. We noticed that steady state T cell count had huge variability. In addition, the T cell model from [Thomas-Vaslin et al., 2008](https://www.jimmunol.org/content/180/4/2240.long) was developed for FVB mice, but C57BL/6J mice was used in experiments in [Boyer et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30905737/). Lymphocyte counts were very different in these 2 types of mouse. For example, lymphocyte count in C57BL/6J mice averages ~4k per microliter of blood; in FVB mice, the count is 1-2k per microliter of blood, with significant gender difference. 

2. Difference in thymus environment. Literature suggested different thymus environments (e.g. epithelialc cells, cytokines) were known to influence T cell development ([Velardi et al., 2021](https://www.nature.com/articles/s41577-020-00457-z)).

## Adapting mouse model for younger mice

Our eventual goal for the mouse model was to predict the ex-vivo gene therapy for ADA-SCID mice. One challenge was that many experiments were carried out on mice that were 14-18 days after birth to avoid potential systematic complications or death.

However, in these very young mice, their naive T cell dynamics were different. They had higher thymic output, probably due to either a more lymphoid-biased differentation (i.e. higher MPP -> CLP rate), or faster proliferation in lymphoid progenitors ([Young et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5068232/), [Baran-Gale et al., 2020](https://elifesciences.org/articles/56221)). 

Sensitivity analysis showed that the steady state peripheral blood T cell concentration was not sensitive to DN/ DP cell proliferation rates (data not shown). Thus, for simplicity, we tune the rate of MPP -> CLP. We use the CD4+ thymocyte as the benchmark. In prior simulation, this number was ~ 300k. In 15-day mice, it was ~780k ([Lee, 2018](https://theses.gla.ac.uk/31002/)). In addition, splenic B cell count was ~25.2M in these young mice ([Blackburn and Kellems, 2005](https://pubmed.ncbi.nlm.nih.gov/15705418/)).

Our parameter scan suggested that increasing MPP -> CLP rate by 2.5-3.5 folds could address the discrepancy in steady state lymphocyte count between a 15-day mouse and 3-month old mice. 

![](img/mouse_15day_ada_scid.png)

By choosing to increase MPP -> CLP rate by 3 folds, lymphocytes at steady state have higher counts. 

| Mouse | naive T per uL blood |  naive B per uL blood |
| ----- | -------------------- | --------------------- |
| young | 5.3k                 | 325                   |
| adult | 1.7k                 | 130                   |


## Adapting mouse model for Adenosine deaminase deficiency Severe combined immunodeficiency (ADA-SCID) mouse

ADA-SCID was caused by mutations in adenosine deaminase (ADA), an enzyme that was reported to be important for dADP, dATP degradation in cells. Lack of ADA resulted in dADP, dATP accumulation and toxicity to cells. ADA was universally expressed, but most active in lymphocytes ([Flinn et al., 2018](https://ojrd.biomedcentral.com/articles/10.1186/s13023-018-0807-5)).

ADA-SCID manifestation in mice included abnormal thymocyte. Double positive (DP) cells was reported to have increased apoptosis, though dynamics of DN cells was reported to be normal ([Whitmore and Gasper, 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00314/full)). In addition, the apoptosis and distribution of peripheral T cells appeared to be normal in ADA-SCID mice ([Apasov et al., 2001](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC209335/)). 
 
ADA-SCID also resulted in impairment in B cell development. The most prominent phenotypes were on activated B cells, while mature naive B cells in the spleen was reported to have a higher death rate ([Whitmore and Gasper, 2016](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00314/full), [Blackburn and Kellems, 2005](https://pubmed.ncbi.nlm.nih.gov/15705418/), [Aldrich et al., 2003](https://www.jimmunol.org/content/jimmunol/171/10/5562.full.pdf), [Prak, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3366416/)).

Thus, for our adaption, we focused on adjusting DP thymocyte and naive splenic B cell death rates. 

To determine DP thymocyte death rate, we used percentage of DN cells amongst thymocytes as the benchmark. We chose this because this appeared to be the most robust readout, while other benchmarks like absolute progenitor/ lymphocytes counts appeared to have lots of variation ([Blackburn and Kellems, 2005](https://pubmed.ncbi.nlm.nih.gov/15705418/), [Blackburn et al., 1998](https://www.jbc.org/article/S0021-9258(17)47085-2/pdf), [Blackburn et al., 2000](https://rupress.org/jem/article/192/2/159/30028/Metabolic-Consequences-of-Adenosine-Deaminase), [Morreti et al., 2021](https://www.nature.com/articles/s41598-021-02572-w)). 

In control mice, DN cell was reported to be less than 5% of total thymocytes; in ADA-SCID mice, this number was ~ 20% ([Morreti et al., 2021](https://www.nature.com/articles/s41598-021-02572-w)). Parameter scan indicated that that DP death rate should be increased to 2 day<sup>-1</sup>. 

![](img/mouse_dpdeathscan.png)

To determine splenic mature B cell death rate, we looked at splenic B cell count. In ADA-SCID mice, splenic B cell count was reported to be ~ 4M. Parameter scan suggested that splenic mature B cell death rate should be increased from 0.032 day<sup>-1</sup> to 0.15 day<sup>-1</sup>. 

![](img/mouse_BMspl_dscan.png)

Steady state lymphocyte in periphral blood was predicted to have a significant drop in ADA-SCID mouse model. This number was more overall consistent with observations in [Blackburn et al., 1998](https://www.jbc.org/article/S0021-9258(17)47085-2/pdf).

| Mouse    | Lymphocyte in periphral blood        |
| -------- | ------------------------------------ |
| Control  | 5.6k per uL                          |
| ADA-SCID |   776 per uL                         |



## Simulation ex-vivo gene therapy on ADA-SCID mice. 

In [Carbonaro et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22833548/), ex-vivo gene therapy was performed on mice that were either on enzyme replacement therapy (ERT) or not. However, we chose to focus on mice on ERT, as study showed that they had almost normal immune system ([Morreti et al., 2021](https://www.nature.com/articles/s41598-021-02572-w)).  In this study, authors isolated bone marrow cells from donor mice, transduced the cells with an average vector copy (VCN) at 0.99 +/- 0.76 (i.e. on average, each transduced cell carried 0.99 transgene) before infusing 5M cells back to the receiver mice.

The receiver mice were conditioned on 200cGy radiation before receiving the cell infusion. 

### Gene therapy simulation

We assumed the conditioning resulted in 30% loss of all dividing/ differentiating progenitors in the bone marrow, and 80% loss in thymus. Fractions of each progenitor cells in transplanted bone marrow cells was estimated based on steady state bone marrow composition. 

Overall, the simulation result was consistent with early data, and timing for zenith and nadir of VCN were overall consistent with observed data. In addition, simulations at steady state were also consistent with observed data. 

![](img/mouse_adascidgtsimul.png)

# Scale HSC differentiation model from mouse to human

Here, we aimed to show how we scale this integrated mouse HSC differentiation model to human. 

Overall, we reverted some cell amplification round, mean residence time, and cell volume back to human. We adjusted CMP amplification time from 16 to 36 to compensate cells being diverged to myeloid lineage. 

![](img/human_full_structure.png)

## Scale myeloid lineage from mouse to human

We increased GMP amplification number from 32 to 128 to compensate for more granulocytes death in human. 

## Scale B cell dynamics from mouse to human

Since it was almost impossible to directly measure the rates of B cell development in human, while the cellular processes are fairly comparable between these two species ([Carsetti et al., 2004](https://pubmed.ncbi.nlm.nih.gov/14962195/), [Ghia et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9785673/)), we kept most of the mouse parameters, but one significant adjustment on the capacity of bone marrow hosting B cells. We scaled it to be 8.4B. This was estimated by linear scaling the capacity based on cell count of bone marrow between mouse and human listed as follows: 

- In mouse, the total bone marrow had ~500M cells. 
    
- In human, we estimated the total cellular count ~ 700B. This was based on the cell density (0.4M per uL, reported in [Mulder, 1974](https://academic.oup.com/ajcp/article-abstract/61/2/199/1770391?redirectedFrom=PDF)) and volume (~ 1.75L, [Hassan and El-Sheemy, 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1079613/#:~:text=An%20adult%20weighing%2070%20kg,haemorrhage%20it%20can%20increase%20sixfold.)).

- In mouse, the B cell progenitor hosting capacity in bone marrow is 6M. 

Thus, the human bone marrow B cell capacity = 6M/500M * 700B = 8.4B. 


In addition, we also increased propreB cell proliferation rate, as human B cell progenitor was observed to proliferate faster in vitro than being transplanted into mouse ([Kraus et al., 2013](https://pubmed.ncbi.nlm.nih.gov/24379121/)).

## Scale T cell dynamics from mouse to human

One of the major difference between mouse and human naive T cell dynamics is their source: in mouse, thymic output is the main source for naive T cell maintenance; in human, naive T cell proliferation in peripheral tissue and blood have a prominant role ([Famili et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5583695/), [Spits, 1994](https://pubmed.ncbi.nlm.nih.gov/8011206/)). Thus, we introduced naive T cell proliferation in periphery blood and tissues while kept most of the rates in T cell development in thymus the same [Thomas-Vaslin et al., 2008](https://www.jimmunol.org/content/180/4/2240.short). In addition, we also incorporated naive T cell entering and leaving lymph nodes. 

Overall the model structure was shown as follows. Major adjustments were

  + Separated lymph nodes from other tissues. We allowed naive T cells to enter peripheral tissues and lymph nodes with different rates. 
  + Introduced naive T cell proliferations in peripheral tissues and blood.
  + Introduced naive T cell death in peripheral tissues. 
  + No naive T cell death in lymph nodes.
  + Same type of naive T cells had the same death rates across organs, but CD4+ and CD8+ naive T cells had different death rates. 
  + CD4+ and CD8+ naive T cells had different rates leaving peripheral tissues. 


In addition, amplification number of CLP after exported to thymus is increased a fold to compensate for human size.

<img src="img/human_naiveT.png" alt="model structure" width="300" >


All the T cell dynamics rates outside thymus were taken from [Sove et al., 2020](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12546). Note one issue was that this paper didn't provide rates for CD4+ naive T cells, so we mostly used CD8+ naive T cell rates, except the rate for naive CD4+ cells leaving peripheral tissues. This value was tuned to make the CD4+:CD8+ ratio in blood ~ 2 ([Benjamin et al., 2018](https://www.sciencedirect.com/science/article/pii/B9780323401395000838)), as showed below. 


![](img/human_cd4_peripheral_exit_rate_scan.png)

# Model validation

Simulated steady state values and their corresponding physiological values were listed as follows. Simulated steady state values were close to literature values.

|    | simulated | observed |
|----|-----------|----------|
| RBC count (per uL) | 3.9M | ~4M ([Pittman, 2011](https://pubmed.ncbi.nlm.nih.gov/21634070/)) |
| Hb (g/dL) | 13.2 | 12-15 ([Pittman, 2011](https://pubmed.ncbi.nlm.nih.gov/21634070/)) |
| Hb in RBC (g/dL) | 366 | ~330 ([Pittman, 2011](https://pubmed.ncbi.nlm.nih.gov/21634070/)) |
| T cell count (per uL) | 662 | 1243 ([Alamooti et al., 2010](https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.20912), note this number is not naive T cell specific) |
| Thymic output (per day) | 32M | 10M - 2700M ([Bains et al., 2009](https://www.jimmunol.org/content/183/7/4329.long)), note the high thymic output happens in infancy |
| B cell count (per uL) | 137 | 101 ([Perez-Andres et al., 2010](https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.b.20547)) |
|Granulocyte count (per uL) | 1.2k | 1-8k ([Yazdi, 2021](https://labs.selfdecode.com/blog/granulocytes/#:~:text=infections%20%5B27%5D.-,Granulocytes%20Normal%20Range,microliter%20(%C2%B5L)%20of%20blood.)) |


We use data from [Ribeil et al., 2017](https://www.nejm.org/doi/full/10.1056/nejmoa1609677) for validation as well. The validation on blood cell/ hemoglobin reconstitution was similar to what we discussed before. 

![](img/human_erythroid_validation.png)


# ADA-SCID model adaption

We use the death rates on DP thymocytes and splenic mature B cells from ADA-SCID mice in this human model. Observed data was obtained from a patient that was not treated with PEG-ADA in [Sauer et al., 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3366410/).

|                    | simulated | observed |
|--------------------|-----------|----------|
| blood T cells (uL) | 149       | 138      |
| blood B cells (uL) | 61        | 0        |

The simulated blood T cell was close to observed number in [Sauer et al., 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3366410/).

# Ex-vivo gene therapy for ADA-SCID patient <a id='ada-scid-infant'></a>

Data was obtained from [Aiuti et al., 2009](https://www.nejm.org/doi/pdf/10.1056/NEJMoa0805817). In our simulation, we focused on patient #1 because this was the only patient had no prior treatment or adversal event during the treatment. 

We assumed this patient's body weight to be ~7.6kg, based on the patient's gender and age. This translated into a total ~65M of infused CD34+ cells. 

The simulation recaptured the reconstitution of granulocytes, B cells and T cells in patient's peripheral blood. 

![](img/human_ada_scid_GT.png)

# Challenges and uncertainties <a id='challenges'></a>

Data availability has been the biggest challenge for model development. This was especially true for progenitor cells between stem cells and fully differentiated ones (e.g. RBC, B cells). In addition, many rates for cell differentiation and proliferation were not reported. Thus, we largely relied on steady state and parameter tuning to provide educated guesses for unreported rates. Here, we list all the parameters that are not directly taken from literature.

- mouse

  + the rate for common lymphoid progenitor (CLP) being exported to thymus
  + division between CLP and B cell and T cell progenitors in bone marrow and in thymus, respectively
  + radioactive-conditioning strength on different cell population
  + breakdown of the transplanted stem cell/ progenitors 
  + difference between young mouse and adult mouse (MPP -> CLP)
  + DP cell (T cell progenitor) and splenic B cell death rate in ADA-SCID mouse 

- human

  + the amplification time of granulocyte-monocyte progenitors (GMP)
  + B cell capacity in bone marrow
  + rate for CD4+ T cell exit peripheral tissues


In addition, there are questions about the suitability of using this model for pediatric population. In our [simulation](#ada-scid-infant), we showed that our model was able to predict immune cell reconstitution on a 6-month old. We were surprised by this, since the HSC proliferation was built on the hypothesis that adult bone marrow can only carry certain amount of stem cells. It is not clear how applicable this hypothesis is in pediatric population. 


# Future work <a id='future-work'></a>

1. Adapting the model to predict ex-vivo gene therapy outcome on transfusion-dependent <img src="https://latex.codecogs.com/svg.image?\beta" title=" " />-thalassemia. This can be achieved by modulating hemoglobin beta chain synthesis and RBC lifespan in the existing model. 

2. Adapting the model to predict ex-vivo gene therapy outcome on pyruvate Kinase Deficiency (PKD). This would require the following adjustment in the model: 

    i. Introduce RET death. 
    
    ii. Decrease RET mean residence time

    iii. Decrease RBC lifespan. 


3. Incorporating additional cell types (e.g. NK cells, platelet) or refine current model (e.g. on granulocyte development). This will allow modelders to utilize information that are commonly mentioned in clinical studies, such as neutrophil engraftment, pletelet engraftment, and VCN in PBMC.

4. Prediction on the minimal effective dose/ conditioning strength for ex-vivo gene therapy that can cross-correct (i.e. an enzyme can be secreted by a genetic-corrected cell and then being utilized by another diseased cell to correct its phenotype; digram shown below). This may be especially relevant when the the ex-vivo gene therapy for lysosomal storage diseases, as lysosomal enzymes are secreted and reuptaken through the same pathway ([Dahms et al., 1989](https://pdf.sciencedirectassets.com/778417/1-s2.0-S0021925818X60724/1-s2.0-S0021925818638256/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjELz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIBlIdRYr%2B9slXPOpOuWrYhvnWUeH2%2Bl4Ll0%2Ff6uED8u0AiEA9BZa1rj2%2FM5whzjyK3rJrqiK7hMoK7FAzGATlnD1sfQqzAQIVBAFGgwwNTkwMDM1NDY4NjUiDOQU73DXLv%2FPCLnWtyqpBDFJ3clcNqlJa26XbTQgevZnO26bX25kQhV0yFsS5ozAbMWaE%2F4v%2BxIqKrfjItCxjiucf4nw7Shgud60YtnXzj6fA6auII0hAmmHfZ4tlhzOy64rDFHWLW6QPXzlWuEZfP5PnHAVXL2VnnK%2Fxk8HSBwbJCsgdOpV8xh5%2FR92wajInOu6pdETAFnrlbntq%2B98HAyCW5o4wkgBhC2R%2Favu8HB7Onp5%2B0i2CT9Voz36JOSOYTJ14oLz6XhmFhtExk3NFBDEcttqagC7NvIBga4UB9V80TiQ0HOlcI7fCT%2FAFzgf%2BcHDHioLnlL3kIeR9nJnhcGaquIiJBYAc1LXjN1LkydpQiTWyIkApPfTutXiBCHcSw6f%2FoHDm7zBKjilGRDDJEQ%2BGE6wKJ5aQ4%2B2S%2BVzK9UzOCt86Li2m1z9CgakjEucegMRZA7%2BQxKMFZa1mQblnmMklI%2FDJ5lXoueTRtj%2BJkowufze2MLf4QFnMUV7JejftI95otFpjjW7uIJ3XPcqHoK3qfyvWqRy4ULaWvL2yQI2BOWJnEV3ymdkLE4OIvYPN0vO5bB%2FlGZ2RA77gF%2FPKrfemjYFxSJ9tzolgWqSy2tvrQXvhZrvuqo7E27t%2FKVdvYK9OJ9v9fCFxEomiZ2Td3Z8nMhB3pWkapp0ufdLiluitM0G2v5nm1qytdACfF3X0pAdff4YQn2AMSa5iN%2BK3xb%2Fii6ZbIs%2FxFpoQUe7yC3NALWKuIqBC%2BUwkuW7nwY6qQFvUSFfuZ5IhP9ZwpPb7mlgz7g6rwcp%2BXykXiiwxOMCTo%2BTobnFP6R692%2BVAwlbi9T2lJZPfbx77nAAXmDPUg8738SjwkGAJ%2BSopcJFTQR8Vj8zftmBYnwB596rfusPNqIMFjDaGmIH331QtDArba9trnFcxMG67jcpC90etrpj0lNTlSXEyQRgdffh14aATU1p6qKCMyrHCoxf3%2FL6117wsCTTsyzpYbpN&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230217T035309Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQJPGCUVU%2F20230217%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=22421b67565b081edf79ae26a843d051ae79ae018788685157644be44f2f1e2d&hash=ac1afcc7ab25ff61ec9a2bba51670f49bda9e47f5be90442c297acf967730f5b&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0021925818638256&tid=spdf-020ae892-9b8f-4c9c-897b-13ee6d62af2c&sid=7ae36e8318120649e588d807e59f7dccdc0bgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f1650580f055a53510600&rr=79ab930ada828c9b&cc=us)). One of the diseases of interest is metachromatic leukodystrophy (MLD) caused by arylsulfatase A (ARSA)-deficiency (ARSA is a lysosomal enzyme). It has been observed that amongst patients this disease receiving ex-vivo gene therapy, having submyeloablative or myeloablative conditioning showed similar level on VCN in PBMC/ enzyme activity ([Fumagalli et al., 2022](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02017-1/fulltext)), further demonstrate the need to investigate the minimal effective dose/ conditioning strength for ex-vivo gene therapy.

![cross-correction diagram](img/Diagram_CrossCorrection.png)

5. Further investigation on using this model for pediatric population. We've discussed this in the uncertainty part of the model. 

6. Sensitivity analysis on the parameters that influences the dynamics of 

# Content <a id='content'></a>

## Files

- ```ExVivo.Rproj``` (R project file; see setup section)

- ```human_lymphoid.Rmd``` (Main script to generate all figures in the human full model with lymphoid and myeloid arms that were scaled from mouse)

- ```mouse_ADA_SCID.Rmd``` (Main script to adapt mouse HSC differentiation model from adult mouse to young mouse and ADA-SCID mouse for gene therapy)

- ```mouse_erythrocytes_Hb.Rmd``` (Main script to generate all figures in developing HSC differentiation model in mouse)

- ```pkgr.yml``` (file to set up R environment; see setup section)

- ```README.md``` (This README file)

- ```SCD_Zheng2021.Rmd``` (Main script to test the implementation of HSC -> RBC differentiation model from [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638) and parameter tuning)

## Folders

- model (where all models files were stored)
    
    - `fullmodel2_scd.cpp`: implementation of [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638) with modifications described [before](#human-scd-implementation) (summarized [here](#human-scd)). 

    - `mouse_erythrocytes_Hb.cpp`: scaled HSC differentiation model in [Zheng et al., 2021](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638) model from human to mouse (summarized [here](#mouse-scd)). 

    - `mouse_erythrocytes_lymphoid_myeloid.cpp`: incorporate myeloid and lymphoid branches to the mouse HSC differentiation model (summarized [here](#mouse-integrated)). 

    - `human_erythrocytes_lymphoid_myeloid.cpp`: HSC differentiation model scaled from mouse to human (summarized [here](#human-integrated)).

    - `human_scd_lymphoid.cpp`: HSC differentiation model (with erythoid, lymphoid, and myeloid branches) scaled from mouse to human with the SCD hemoglobin dynamics (summarized [here](#human-integrated); note this model tracks hemoglobin slightly differently to adjust for trangene on hemoglobin on this therapy).

- img (where output image files are saved)

- data (source data from publications. See readme.txt in the folder for more details.)

- julia (julia implementation of the integrated model in human. See README inside the folder for more details)

# Setup<a id='r-setup-section'></a>

## Prerequisites

R version 4.1.3 and the following R packages: 
- mrgsolve
- tidyverse
- gridExtra
- grid
- ggplot2
- mrggsave

## package management

1. Install pkgr software following the instructions [here](https://github.com/metrumresearchgroup/pkgr). 
2. Open the R project `Public.Rproj`. This allows you to work from within a self-contained project environment.
3. Install packages by typing in terminal: `pkgr install`. This command will look for the file `pkgr.yml` and install the packages listed. The specific package versions are imported from https://mpn.metworx.com/docs/.

## mrgsolve installation
For detailed instructions on mrgsolve installation and important dependencies, follow [this link](https://github.com/metrumresearchgroup/mrgsolve). 
