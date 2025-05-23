---
title: Conclusions
date: "2022-09-03"
output:
  pdf_document:
    highlight: tango
numbersections: true
bibliography: /home/griessbaum/Dropbox/UCSB/library/library.bib
link-citations: yes
linkcolor: blue
header-includes:
    - \usepackage{hyperref}         
    - \usepackage{siunitx}
    - \DeclareSIUnit\month{month}
    - \usepackage{lmodern}
    - \usepackage[acronym, nomain]{glossaries}    
    - \include{tex/acronyms.tex}
    - \usepackage[outputdir=chapters_standalone]{minted}
    - \makeglossaries
---

File-centric data analysis is a paradigm in which files are the smallest unit of data. While file-centricity simplifies the task of archiving and distributing data, it pushes the burden of extracting, transforming, and loading (ETL) data before performing any data analysis to the data users.
Data-centric data analysis, on the other hand, is a paradigm in which the smallest unit of data are individual observations (e.g., instances/objects in some form of a schema), and data thus are stored in some form of a 
database and accessed by some form of query. 

I addressed two questions that need to be solved to allow data repositories and data users to move toward data-centricity:

In chapter \ref{Chapter_OCCUR}, I described (+OCCUR); a system providing identity and citations to data in a data-centric world. 
Identifying data in a data-centric world is different from identifying data in a file-centric world: In a file-centric world, data is accessed through, e.g., file paths or (+^URL), and files themselves can be understood as the identities of the data they contain. In a data-centric world, the notion of files does not exist. Data is instead accessed through queries; therefore, we must be able to identify the result sets of queries.
(+OCCUR) allows us to identify data that is queried through (+OPeNDAP). However, it is intended to be a reference implementation that can be adapted to any other data distribution system 
in a data-centric world. That is, any data accessed through some form of a query can be identified by a system similar to (+OCCUR).

In chapter \ref{Chapter_Software}, I provide a solution to harmonize spatial data. Harmonizing data means ensuring that things that are the same are referred to as the same. Only if data is harmonized can we associate data from different datasets and perform data analysis across multiple datasets. Harmonizing data, therefore, is a requirement to void the necessity of ETL and fully leverage the benefits of data-centricity. In the context of spatial data, harmonization means that there needs to be one unified method to express location. The concept of (+STARE) provides such a unified method, and the (+STARE) software collection implements the (+STARE) concept. Any spatial object can be represented within (+STARE). Further, evaluating spatial relations between spatial objects represented in (+STARE) representation is cheap, allowing one to associate different datasets by their location. While the world may remain in a transitional phase between file-centricity and data-centricity for some time, the (+STARE) software collection provides a bride toward data-centric spatial data analysis. With its rich capabilities to load and convert conventional spatial representation formats, the (+STARE) software collection provides a method to harmonize any spatial data, e.g., at the beginning of a data processing pipeline, and thus allow for a data-centric workflow. 

In chapter \ref{chapter_3}, I provide an example of how scientific data analysis of remotely sensed data can be performed in a data-centric world and highlight the benefits of performing spatial data analysis on spatially harmonized data with the (+STARE) software collection. Previous methods to harmonize spatial data have been based on location discretization and data sampling using regular grids. Those methods fail us when working with data collected at varying spatial resolutions since they require us to re-grid and re-sample data, both of which likely entail suboptimal compromises. Further, the location discretization reduces data's spatial fidelity, resulting in noise in derived products. By working with data that has been harmonized through (+STARE) instead, I avoided both of these problems. By avoiding spatial discretization, I increased the accuracy of an algorithm for fractional snow cover estimations and reduced the noise in the time series of those fractional snow cover estimations. By having all input data harmonized through STARE, I could effortlessly spatially associate grid cells, regions of interest, and individual observations from WorldView Legion, (+MODIS), and (+VIIRS). Other algorithms and data analysis efforts likely will benefit from working with data that has been spatially aligned through (+STARE), both in terms of simplifying the ETL process and in terms of using the full spatial fidelity of observations.
