---
title: Towards the twilight of file-centricity
author: Niklas Griessbaum
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
    - \usepackage[acronym]{glossaries}    
    - \makeglossaries
    - \include{context/acronyms.tex}
    - \usepackage{tipa}
include-before:
    - \printglossary
---

\clearpage

# Motivation
Environmental informatics is the application of information technology to environmental sciences [@Frew2012]. As such, it addresses the information infrastructure that environmental scientists leverage to obtain knowledge from environmental data. 

The appearance of cloud computing, characterized by scalability on one side, 
and by abstraction and stereotyping of interactions through services and (+^API) on the other side [@Foster2017] provides for opportunities in environmental informatics. It opens up questions about what the workflow of environmental scientists in the 21st century should look like to fully exploit these opportunities.

My thesis is that the inherent heterogeneity of the flow from data to knowledge in the environmental sciences causes bottlenecks that can be unblocked by moving processing paradigms away from file-centricity and towards data-centricity. In my dissertation, I address the "twilight of file-centricity"[^dateidaemmerung] and technologies required to transition from file-centricity to data-centricity.

[^dateidaemmerung]: German: Dateidämmerung \textipa{[\textsubring{d}a'ta\textsubcircum{i}'dEm@\;RUN]}} 

Files package (i.e. chunk or aggregate) data into logical units and provide intelligible identities to their contents: their filenames. For example, MODster [@Frew2005; @Frew2002] used the conventions of (+MODIS) granule filenames as data identity to managed distributed data. Files further conceal the structure of data they are holding, which allows repositories to preserve and distribute data structure-agnostically. While this generalizes and simplifies the task for data repositories, it pushes the responsibility of acknowledging the data's structure to the users, who have to (+ETL) data prior to extracting knowledge [@Rilee2016; @Szalay2009].

The problem in this approach materializes in two bottlenecks: data movement and data alignment.

## Data Movement:
A system that is agnostic of the structure of the data it is holding is incapable of performing computations on the data. It can merely act as a point of preservation and distribution. Data therefore has to be moved (more accurately: copied) to the point of computation (Gray's 3rd Law, [@Szalay2009]).

Data movement however is undesired since it results in uncoordinated and unstructured duplication of data and therefore in storage waste.
Further, we are facing an increasing disparity between network speeds and compute power. User bandwidth speeds have been following Nielsen's Law [@Nielsen1998] and grew annually by \SI{50}{\percent} over the last 36 years while compute power has been following Moore's law [@Moore1975] and grew annually by \SI{60}{\percent} for the last 40 years[^moore],  making it less and less attractive and ultimately infeasible to move data to the point of computation as data volumes grow [@Hey2009]. While researchers may choose to copy gigabytes worth of data for their analysis, copying petabytes will not be an option within the near future [@Szalay2006].

[^moore]: [https://ourworldindata.org/grapher/transistors-per-microprocessor](https://ourworldindata.org/grapher/transistors-per-microprocessor)

The predefined package size of files makes matters worse: If the package size does not exactly equal the area of interest for an analysis[^file_package], a transfer overhead is incurred [@Gray2002].

[^file_package]: A file might for example contain bands, areas, or time periods not needed for a given analysis.

## Data Alignment
Files allow repositories to ignore data structures; they therefore free repositories from the responsibility of data alignment. It becomes the task of the data user to align data during the (+ETL) process if various datasets are to be integrated. Data alignment provides a common method to express (e.g., spatiotemporal) coincidence throughout all datasets. It can be achieved through (re-)projection, aggregation and gridding, and/or indexing. Because methods for alignment are not provided from a central point (i.e. the repository), every user has to manually align data, resulting in redundant effort.

In order to break the stated bottlenecks, the file-centric approach has to be replaced by a system that:

- Is capable of providing identity to data independently of packaging.
- Does not require manual data alignment.
- Voids the necessity to transfer data by being aware of the data structure and therefore being capable of performing computations at the place of storage.

\newpage

# Summary

In my dissertation, I address the following challenges to the twilight of file-centricity.

## Data Identity and Data Citations
By abandoning files as the package of data, a natural identity of data is lost. Identity, however, is needed to refer to and reason about data. Regardless of previous subsetting and processing, data needs to be identifiable and citable. I therefore developed a service providing identity and simplifying the creation of citations to data stored in, and served through, online repositories.

Citations help to make data (+FAIR) [@Wilkinson2016]. In more abstract terms, data citations provide identity to data, which allows referencing and de-referencing. Provision of identity is a key challenge in the twilight of file-centric workflows, since a natural addressable identity of data is lost as soon as files as a package of data are abandoned.

Data citations differ from citations of printed material in that the cited content (i.e. the data) may evolve over time and in that meta-information such as authorship or the provenance may vary within a continuous dataset [@Buneman2016]. Further, since generally speaking infinite ways to subset data are possible, data citations cannot be statically generated. They have to be machine-actionable, both in terms of dynamic creation (as a function of time and subsetting parameters) and in terms of resolving data citations back to the cited material.

Technologies such as the (+WCS) and the (+OPeNDAP) [@Gallagher2005] play a key role in the twilight of file-centricity. Their ability to seamlessly provide access to data rather than to files provides an ideal starting point for theoretical and practical excursions on how to address identity and citations in a data-centric workflow. 

With the development of the web service (+OCCUR)[^occur], I am exploring an approach for assigning identity and citations to dynamic data. (+OCCUR) is a web service that allows users to assign and store identities for data retrieved from an (+OPeNDAP) query. (+OCCUR) creates identifiers for identities which can later be resolved through (+OCCUR), whereby (+OCCUR) will verify that the data has not changed since the assignment of the identity. (+OCCUR) further brokers identities by ensuring that identical data is sharing the same identity.

[^occur]: [https://github.com/NiklasPhabian/OCCUR](https://github.com/NiklasPhabian/OCCUR), [http://occur.duckdns.org](http://occur.duckdns.org) 

(+OCCUR) additionally allows users to generate formatted citation snippets for both (+OCCUR) identities and for any (+OPeNDAP) query. It is the (+OPeNDAP) counterpart to the \textit{DOI Citation Formatter[^crosscite], which allows the creation of citation snippets from (+DOI). 

[^crosscite]: [www.crosscite.org](www.crosscite.org)

The development of (+OCCUR) was supported by the (+ESIP) federation as a 2018 (+ESIP) lab project. The work on (+OCCUR) has been presented at the 2017 Bren PhD Symposium and the 2018 (+ESIP) summer meeting.

## Data alignment
To make ETL processes obsolent, potentially heterogeneous data have to be co-aligned. Co-alignment means that a common concept to address coincidence throughout all datasets exists. It can be achieved through a common indexing schema. A side-effect of co-alignment is that it can be exploited to improve physical data placement: The index can be used to store coinciding data in physical proximity. This will, for a lot of use-cases, allow processing data in parallel in a shared-nothing architecture [@Kuo2017]. The prevalent dimensions for which data has to be aligned in environmental science are time and space.

I am addressing data co-alignment with the development of the STARE software collection. The (+STARE) is a global indexing schema that is built upon the quadtree (+HTM) [@Kuo2017] and provides a common concept to address spatiotemporal coincidence for environmental data.


The STARE software collection is a first step towards true data centricity. STARE is a central technology allowing all datasets required for a given analysis to be co-located at the place of computation. STARE-indexed geospatial data 

The co-located datasets further need to be co-aligned to allow interoperability [@Kuo2017; @Rilee2016] and make (+ETL) processes unnecessary. In practical terms, this means storing the data in a database.
 

I have presented the findings of an early stage of EarthDB 2.0 at the 2017 (+AGU) Fall meeting. I want to continue the development of EarthDB 2.0 to create a system capable of solving common problems in environmental science. In particular:

\paragraph{STARE in Python (pystare)}
Though the existence of STARE may be completely opaque to a user of EarthDB 2.0, exposing STARE functionality to a high-level programming language will be required for educational and debugging purposes. I therefore want to enable common Python geometry representations such as \textit{shapely} and \textit{geopandas} to be convertible to and from STARE representations.

\newpage


## Use Cases
One of the five laws postulated by Jim Gray is that a database designed for a given discipline has to be capable to answer the 20 key questions a scientists may have [@Hey2009; @Szalay2009]. Following this spirit, I evaluate the previously proposed system with two remote sensing use cases: time series of night lights, and multi-sensor snow mapping.  

In order to verify the usability of the previously proposed solution of EarthDB 2.0 in environmental science, I am going to use EarthDB 2.0 to solve environmental science domain use cases. Besides solving a domain problem itself, the work on these uses cases will enable me to exhibit typical interactions between environmental scientists and EarthDB 2.0.

I am intending to carry out two distinct use cases, the first one to verify general and simple spatiotemporal functionality, and the second one to demonstrate that complex environmental uses cases can be addressed with EarthDB 2.0.

### Night lights
The (+VIIRS) onboard the Suomi (+NPP) has a (+DNB) that is sensitive to visible and near-infrared wavelengths, which enables it to observe nighttime lights on Earth at a significantly higher spatial and temporal resolution than its predecessors such as the (+DMSP)-(+OLS). 

Prof. Mark Buntaine proposed to use timeseries of averaged nighttime light intensity during/after natural disasters over a given administrative area as a predictor of the area's resilience to disasters. 
As the (+VIIRS) (+DNB) product is only available as swath data (and hence each pixel is located in an irregular grid), a conventional (+GIS), such as PostGIS, would require billions of point-in-polygon tests in order to associate individual DNB pixels with an administrative area. Even with R-tree indexing of the administrative areas, an on-the-fly workflow is impossible even for a relatively small spatial extent (e.g., Puerto Rico.)

In order to demonstrate EarthDB 2.0's capability, I want to import \textit{all} (+VIIRS) (+DNB) data for the continent of Africa and enable users to extract timeseries of averaged nighttime light intensities for arbitrary polygons (e.g., each level 2 administrative area) on-the-fly. 

I am intending to demonstrate my first findings of this approach during the 2020 ESIP winter meeting and refined findings during the 2020 ACM SigSpatial meeting.


### Multi sensor snow mapping
I will implement a snow mapping system that uses low-level/swath data from multiple sensors, for two reasons. Firstly, it allows for cross verification of results. This is a requirement for the development and testing of any novel algorithm. Secondly, the use of multiple sensors increases the aggregated revisit rates, therefore increasing the temporal and/or spatial resolution of snowmaps in the sense of (+HRPP), which are designed to provide the best precipitation estimate at any given time using data from multiple satellites [@Lettenmaier2015].

The use of lower level products will avoid artifacts caused by resampling and further increase the spatial and temporal resolution by circumventing the pre-defined gridding of higher level products. High level products use a referencing matrix that provides information on coordinates of an image corner, pixel spacing, and rotation and thereby allows calculation of the coordinates of any point in the image. However, translating imagery acquired from satellites into such a coordinate system requires resampling of the measured radiances, thereby risking introduction of artifacts into the data.}. Additionally, it will allow for the integration of products that are only available at lower levels, such as the thermal bands from the (+MODIS) calibrated radiance product (MOD02KM / MYD02KM).

Multi-sensor sub-pixel snowmapping is an interesting use case to test EarthDB 2.0 for two reasons. Firstly, (+STARE) is intended to facilitate the integration of inhomogeneous data from various sensors at different spatial and temporal resolutions. Secondly, (+STARE) is intended to facilitate the integration of lower level swath products, which are otherwise cumbersome to work with. Both of these assumptions can be tested in this use case.

\paragraph{Use of swath data:}
(+MODSCAG), as described by [@Painter2009], is an algorithm based on linear endmember unmixing to retrieve sub-pixel snow cover data as well as grain size and albedo estimates from MOD09GA (level (+L2G). L2G is the "Level 2G" format, which is geolocated, and gridded into a map projection.) data. 
Using gridded data allows ignoring issues connected to the viewing geometry of MODIS, causing a) variations in pixel sizes in scan directions and b) anisotropic reflectances at shallow viewing angles. I propose to relax these simplification by porting (+MODSCAG) to use (+MODIS) swath data (MOD09) within EarthDB 2.0.

\clearpage

# References {-}