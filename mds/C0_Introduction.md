---
title: Towards the twilight of file-centricity in environmental sciences
author: Niklas Griessbaum
date: "2022-09-03"
output:
  pdf_document:
    highlight: tango
numbersections: true
bibliography: /home/griessbaum/Dropbox/UCSB/library/library.bib
suppress-bibliography: false
link-citations: yes
linkcolor: blue
header-includes: 
    - \usepackage{hyperref}         
    - \usepackage{siunitx}
    - \DeclareSIUnit\month{month}
    - \usepackage{lmodern}
    - \usepackage[acronym]{glossaries}
    - \input{tex/acronyms.tex}
    - \makeglossaries
    - \usepackage{tipa}
include-before:
    - \printglossary
---

\clearpage

# Motivation {#c0_introduction}
Environmental informatics is the application of information technology to environmental sciences [@Frew2012]. As such, it addresses the information infrastructure that environmental scientists leverage to obtain knowledge from environmental data. 

Environmental data is traditionally collected, archived, and distributed in computer files. Alike their real-world counterpart, computer files are containers holding content (i.e., the data) and have intelligible labels (i.e., filenames). Since files are mere containers, there is no prescribed structure for their content. Reading files, therefore, requires contextual knowledge (aka metadata) to interpret what the content means.

Packaging data into files generalizes the tasks of archiving and distributing data. The discretized nature of files makes it easy to reason about their identity: We may state where a file is located, what its name is, evaluate its size, and even compute checksums. Further, files can be archived and distributed without considering the structure of the content or how the content is to be interpreted. While this simplifies the task for data repositories, it pushes the responsibility of acknowledging the data's structure to the users, who have to (+ETL) data prior to extracting knowledge [@Rilee2016; @Szalay2009].

> **My thesis is that the content-structure-agnostic nature of files causes unnecessary bottlenecks in the flow from data to knowledge in environmental sciences. Unblocking those bottlenecks requires moving data processing paradigms away from file-centricity and towards data-centricity. In my dissertation, I address the "twilight of file-centricity"[^dateidaemmerung] and technologies required to transition from file-centricity to data-centricity.**

[^dateidaemmerung]: German: Dateidämmerung \textipa{[\textsubring{d}a'ta\textsubcircum{i}'dEm@\;RUN]}

# File-Centricity 

The issue with file centricity is that it incurs redundant work and prohibits data processing at the point of storage: A system agnostic of the structure of the data it holds is incapable of performing computations on the data. It can merely act as a point of preservation and distribution. Data, therefore, has to be moved (more accurately: copied) to the point of computation (Gray's 3rd Law, [@Szalay2009]). However, data movement is undesired since it results in uncoordinated and unstructured data duplication and storage waste[^datamove]. But it is not merely the duplication of the data that is problematic. Moreover, every user that copies data must redundantly acknowledge the data's structure during (+ETL).

[^datamove]: Further, we face an increasing disparity between network speeds and compute power. User bandwidth speeds have been following Nielsen's Law [@Nielsen1998] and grew annually by \SI{50}{\percent} over the last 36 years, while compute power has been following Moore's Law [@Moore1975] and [grew annually by \SI{60}{\percent} for the last 40 years](https://ourworldindata.org/grapher/transistors-per-microprocessor), making it less and less attractive and ultimately infeasible to move data to the point of computation as data volumes grow [@Hey2009]. While researchers may choose to copy gigabytes worth of data for their analysis, copying petabytes will not be an option within the near future [@Szalay2006]. The predefined package size of files makes matters worse: If the package size does not exactly equal the area of interest for an analysis (A file might, for example, contain bands, areas, or periods not needed for a given analysis.), a transfer overhead is incurred [@Gray2002].

During (+ETL), users align data from various sources. Aligning data means harmonizing attributes and dimensions across datasets, allowing us to associate and compare data. In other words, alignment means that a common concept to address coincidence throughout all datasets exists. In the environmental sciences, characterized by a prevalence of spatiotemporally resolved data from observations and models, it is often of interest to associate spatiotemporally coincident data and thus to spatiotemporally align data. 

However, spatiotemporal alignment is cumbersome and challenging. There is a multitude of concepts, file formats, and referencing scheme in which spatiotemporal data is stored and expressed: To name a few, locations may be conceptualized as continuous fields or as discrete features, expressed through affine transformations or as lists of coordinates, stored as projected image files or in relational databases. Time, and more so time duration, may be expressed by many calendrical formats or as offsets of epochs with or without consideration of leap seconds. Consequently, spatiotemporal alignment is typically an error-prone tailormade process involving a multitude of compromises and much work.

# Towards Data-Centricity 

Contrary to file-centricity, data-centricity requires data to be co-aligned. Data alignment voids the necessity of (or at least simplifies) (+ETL) and makes computation at the point of storage possible[^compute_storage]. In practical terms, this means storing data in some kind of database[^database]. 

[^compute_storage]: The ultimate goal of data-centricity is voiding the need (or at least reducing) data movement. Data alignment is hereby crucial: Having data aligned allows for improved data sharding/placement. I.e., spatiotemporally coincidental data can be stored in physical proximity in, e.g., shared-nothing architectures.

[^database]: “A structured set of data held in computer storage and typically accessed or manipulated by means of specialized software.” [@databaseOED].

I am addressing three questions arising in the file-twilight of environmental sciences:

#### How do we handle data identity in a data-centric world? {-}
Citations help to make data (+FAIR) [@Wilkinson2016]. In less abstract terms, data citations provide identity to data, which allows referencing and de-referencing. Provision of identity is a critical challenge in the twilight of file-centric workflows since a natural addressable identity of data is lost as soon as files as a package of data are abandoned.

Technologies such as the (+WCS) and the (+OPeNDAP) [@Gallagher2005] play a vital role in the twilight of file-centricity. Their ability to seamlessly provide access to data rather than to files provides an ideal starting point for theoretical and practical excursions on how to address identity and citations in a data-centric workflow. 

With the development of the web service (+OCCUR)[^occur] (chapter \ref{Chapter_OCCUR}), I am exploring an approach for assigning identity and citations to dynamic data. (+OCCUR) is a web service that allows users to assign and store identities for data retrieved from (+OPeNDAP) queries. (+OCCUR) creates and stores identifiers for identities which can later be resolved through (+OCCUR), whereby (+OCCUR) will verify that the data has not changed since the identity assignment. (+OCCUR) further brokers identities by ensuring that identical data shares the same identity.

[^modster]: For example, MODster [@Frew2005; @Frew2002] used the conventions of (+MODIS) granule filenames as data identity to manage distributed data. 

[^occur]:  [http://occur.duckdns.org](http://occur.duckdns.org). [https://github.com/NiklasPhabian/occur](https://github.com/NiklasPhabian/occur)

#### How do we spatiotemporally align data to enable data-centric environmental science? {-}
Time and space are the most prevalent (and, simultaneously, most challenging) dimensions that must be aligned in the environmental sciences. We, therefore, require a universal method to express space and time. 

I am addressing data alignment with the development of the (+STARE) software collection (chapter \ref{Chapter_Software}). (+STARE) is a spatiotemporal referencing and indexing schema that is built upon a (+HTM) quadtree [@Kuo2017] and provides a common concept to evaluate  spatiotemporal coincidence for environmental data. The (+STARE) software collection is a first step towards true data-centricity since it allows all datasets required for a given spatiotemporal data analysis to be stored in the same spatiotemporal (database) schema. The (+STARE) software collection contains software to convert conventional files containing spatiotemporal data into the (+STARE) schema and provides a variety of data-centric storage backends.

#### How do we perform environmental science in a data-centric world? {}

Jim Gray's 4th Law [@Hey2009; @Szalay2009] postulates that a data engineering challenge should be approached by determining the 20 most important questions a researcher may want a given data system to answer. Following this spirit, I solved a set of science use cases to drive the development of the (+STARE) software collection. Chapter \ref{Chapter_Software} contains two smaller undertakings: In section \ref{night_lights}, I generate time series of night lights from (+VIIRS) data for a set of administrative areas to determine the characteristics of night light intensity drop and recovery during and after natural disasters. In section \ref{movingobjects}, I track precipitation events and extract spatiotemporal incident data from various sensors. Finally, chapter \ref{chapter_3} provides an approach to improve the accuracy of (+fSCA) retrievals. I here exploit the relatively high spatial accuracy of (+MODIS) geolocations and use (+STARE) to align irregularly spaced observations.

All three use cases demonstrate how (+STARE) allows for integrating inhomogeneous data from various sensors at different spatial and temporal resolutions by providing a harmonized schema for time and space and thus allowing for data-centric workflows.


\clearpage

# References {-}
