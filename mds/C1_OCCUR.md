---
title: OCCUR - An automated data citation system for OPeNDAP resources of remote sensing data
author: Niklas Griessbaum, James Frew
date: "2022-09-03"
output:
  pdf_document:	
    highlight: tango
numbersections: true
graphics: true
bibliography: /home/griessbaum/Dropbox/UCSB/library/library.bib
link-citations: yes
linkcolor: blue
header-includes:
    - \usepackage{hyperref}         
    - \usepackage{siunitx}
    - \DeclareSIUnit\month{month}
    - \usepackage{lmodern}
    - \usepackage[acronym, nomain]{glossaries}    
    - \input{tex/acronyms.tex}
    - \makeglossaries	
	- \usepackage[outputdir=chapters_standalone]{minted}
---

\clearpage


# Abstract {- #c1_abstract}
(+OCCUR) is a web service that creates identifiers and citations for data served by (+OPeNDAP) servers. As a partial implementation of the (+RDA) (+WGDC) guidelines, it addresses the need to identify arbitrary subsets of revisable datasets. (+OCCUR) creates identifiers from a combination of an (+OPeNDAP) query and timestamp and saves a hash of the query's result set. (+OCCUR) can then dereference these identifiers to access the data via (+OPeNDAP) and generate human-readable citation snippets. When accessing data via an identifier, (+OCCUR) compares the saved hash with the hash of the retrieved data to determine whether the data has changed since it was cited. (+OCCUR) uses CiteProc to generate citation snippets from the identifier's (+OPeNDAP) query, timestamp, dataset-level metadata provided by the (+OPeNDAP) server, and optionally the result set hash.

\clearpage
\glsresetall

# Introduction {#c1_introduction} 
## Why Data citations
[@Hey2009] coin the term *fourth paradigm* as "using computers to gain understanding from data created and stored in our electronic data stores." The fourth paradigm adds the exploration of data collected from instruments and simulations to the traditional empirical, theoretical, and computational approaches to scientific research. In this context, data collection and assembly are themselves significant research activities [@Frew2012].

A crucial step into the fourth paradigm is acknowledging data as first-class research products. As such, data must be persistently available, documented, citable, reusable, and possibly peer-reviewed [@Callaghan2012; @Kratz2014] - a process summarized as making data "(+FAIR)" [@Wilkinson2016]. Data citation is one of the required building blocks to achieve this goal, and widespread adoption of data citations is expected to benefit the progress of science [@CODATA2013; @JDDCP2014].

However, there is no consensus on what data publication means [@Kratz2014] nor on how data citation mechanisms are to be implemented [@Costello2009]. The lack of data citation standards was criticized more than a decade ago by [@AltKin07]. Years later, [@Altman2015] and [@Tenopir2011] found that even though required by publishers, researchers still too often do not make data publicly available nor cite the data consistently. There are both cultural and technical reasons for this.

[@Lawrence2011] find that traditionally only conclusions are valued; little attention is given to the fitness of the data for *re-* interpretation. This reduces the motivation for data production and publishing. Even more so, [@Tenopir2011] stresses that researchers may be motivated to purposely withhold data to retain their own ability to publish findings.

On the technical side, robustness, openness, and uniformity in data publication still need to be improved [@Starr2015; @Koltay2016]. Cost, not so much for the storage but for curation efforts, is another reason preventing data publication [@Gray2002]. [@Tenopir2011] states that a significant reason for data withholding is the effort required to publish data. Additionally, [@Belter2014] finds that data citation practices are inconsistent even when used. [@Assante2016] illustrates citation practices ranging from exporting a formatted citation string (or a generic format such as (+RIS) or BibTex[^bibtex]) to embedding links to the data to sharing data on social media.

[^bibtex]: [http://www.bibtex.org/](http://www.bibtex.org/)

To address the need for more standards and uniformity, The (+RDA) (+WGDC) released 14 recommendations to enable automized and machine-actionable identification and citation of evolving and subsettable datasets [@Rauber2015a; @Rauber2015]. [@Rauber2021] describes several (partial) implementations of those recommendations.

[@Silvello2017] provides an exhaustive review of data citations' current state in terms of motivations and implementations. Based on a meta-study, the author identified six motivations for data citations: Attribution, connection, discovery, sharing, impact, and reproducibility. We simplify these to identity, attribution, and access:

### Identity
Citations provide an identity to data, enabling referencing and reasoning about data [@Bandrowski2016] even in its absence (e.g., no longer extant or inaccessible behind a paywall). Identifying data also allows for evaluating its usage, relevance, and impact [@Honor2016].

### Attribution
Citations attribute data to authors, allowing them to take professional credit for it and providing accountability to the sponsors of the data's collection/creation and publication. This provides an incentive for sharing [@Niemeyer2016; @Callaghan2012; @Kratz2014].

### Access
A citation provides information on how to retrieve the cited material (e.g., the journal, year, and pages). Persistent access to data is essential to enable reusability and reproducibility [@Starr2015]. 

<!-- 
Data citations may provide an additional side effect in improving data accessibility: Evidence suggests that well-curated data will be cited favorably [@Belter2014]. Therefore, citations may motivate researchers to curate their data, thus making it more accessible. -->

## What is data citation?
Citations provide identity, Attribution, and access mechanisms to cited material. Data citations differ from citations of printed material in that the cited content (i.e., the data) may evolve and in that meta-information such as authorship or provenance may vary within a continuous dataset [@Buneman2016]. Further, data citations cannot be statically generated for subsettable data unless the number of possible subsets is trivially small. This is specifically true when data, rather than files, are accessed through, e.g., (+^API)[^webapi]. Data citations, therefore, have to be machine-actionable, both in terms of dynamic creation (as a function of time and subsetting parameters) and in terms of resolving citations to the cited material [@Assante2016; @Altman2015; @Buneman2016].

[^webapi]: E.g. (+WMS), (+WCS), or (+OPeNDAP) [@Gallagher2005]

Data citations often use actionable (+^PID) such as (+^DOI). However, actionable (+^PID) blurs the distinction between identity and access [@ESIP2012a]. In this context, [@Buneman2010] emphasizes that (+^DOI) should be considered a part of, but not a substitute for, data citations. Identity and access remain two distinct facets of a citation, and there is utility in data identity regardless of whether or not the data can be accessed or even still exists. [@Parsons2013a] criticize another aspect of (+DOI) use in data citations: (+^DOI) are misunderstood to provide imprimaturs and persistence. However, a (+DOI) cannot provide persistence and should solely be understood as a locator and identifier, which is required long before an imprimatur can be issued.

In the following, we address questions related to data citation:

1. How are datasets and their subsets identified?
2. How is fixity assured?
3. How are revisable datasets handled?
4. How do citations facilitate access to data?
5. How are human-readable citation snippets/strings generated?

Answers to these questions vary widely depending on the scientific domain, as well as particular dataset characteristics such as complexity (tables, arrays, graphs), volume, update frequency, and the repository's services, specifically regarding subsetting.

### Identity
A unique identity can be represented by any arbitrary unique string [@CODATA2013]. In some contexts, already established identifiers such as filenames [@Buneman2016] or accession numbers [@Bandrowski2016] may serve this purpose. In practical terms, [@AltKin07] suggests that the identity should double as a handle to the data by associating it with a naming resolution service, e.g., through the use of, e.g., a (+DOI), (+LSID), (+URN), (+HNR)[^3], or (+URL).

[^3]: [http://www.handle.net/](http://www.handle.net/)

Contrary to traditional publications, data may be queried to produce a potentially unlimited number of subsets from a single source [@Davidson2017; @CODATA2013]. It is, therefore, necessary to reference a dataset and every possible subset of a dataset. [@AltKin07] coin the term "deep citation" to describe the ability to reference subsets of data. Further, data may evolve, which opens the discussion of how to identify the varying states of a dataset [@Huber2015].

This further prompts the question of at which granularity a unique identity should be assigned. [@Buneman2010], therefore, introduce the concept of a "citable unit": An object of interest such as a fact stated in a scientific paper or a subset within a dataset.

<!---The authors argue against creating identifiers for every object of interest but rather to choose the granularity of the citable unit wisely. A citation of an object can then be created by appending the identity of the citable unit with information about the object's location within the citable unit. For example, relational databases could define citable units through views [@Buneman2016]. --->

<!-- [@Rauber2015] and [@Proll2013], as well as the RDA WGDC [@Rauber2015a], suggests identifying subsets by storing (normalized) queries and associating them with a PID. <!--  Since methods and syntax for subsetting depend on the individual repository technology, [@AltKin07] suggests citing the entire dataset and appending the citation with a textual description of how the dataset was subsetted. The authors further suggest that if significant pre-processing is done on the dataset, it should be stored and cited with its own identity individually. --> 

### Fixity
Datasets may change unintentionally through malfunctions, malicious manipulation, or intentionally due to (+CUD) operations. Data fixity is the property of data to remain unchanged, and fixity checking verifies that data has not changed. If fixity checking is included in a data citation system, it can verify if a data object is identical to the referenced data.

[@AltKin07; @Rauber2015; @Crosas2011] suggest including (+^UNF) into data citations to allow fixity checking. Though mentioned in most discussions of data citation systems (e.g. [@Buneman2016; @Davidson2017]), few data citation system implementations have addressed fixity so far.

<!-- In this context, [@Klump2016] reports that the STD-DOI[^4] metadata schema rejected the inclusion of UNFs as they would be dependent on the form of representation (e.g., file formats or character encoding). [^4]: <abbr>DFG</abbr> Project Publication and Citation of Scientific Primary Data-->

### Revisability
Datasets may evolve through updating, appending, or deletion [@Klump2016]. We will refer to these datasets as *revisable data* to distinguish them from datasets changing due to errors or malicious manipulation.

Literature frequently intermixes the term "fixity" and the ability to cite revisable datasets. A revisable dataset is anticipated and intended to change its *state* over time. A citation system consequently has to be able to distinguish between the states of a revisable dataset [@Rauber2015; @Klump2016]. However, to achieve this, merely the abstract state a citation is referencing has to remain fixed (i.e., there cannot be an ambiguity of the referenced state). This is true independently of the ability to dereference a citation to the referenced state (i.e., the actual state being fixed). Identifying and referencing a dataset's ephemeral state is required for data citation. However, the ability to persistently retrieve this state is a data publication, not a data citation challenge. It is up to the publisher to choose an apt level of zeal:

- **Pessimistic:** data is assumed to be ephemeral; consequently, citations cannot ever be dereferenced.
- **Optimistic:** data is assumed to be fixed. Citations always dereference to the current state of the data.
- **Opportunistic:** data is assumed to remain fixed for some time. Citations can be dereferenced only until the data changes.
- **Pedantic:** Every data state is saved; consequently, citations can always be dereferenced to the referenced state.

<!-- [@Klump2016] make a hard distinguishing between growing and updated datasets. Retrieving any state of a growing dataset can be implemented through time ranges, given that records are timestamped.-->

Versioning can be used to identify states of revisable data. We acknowledge that the term *version* is often used synonymously with the term *state* (of a dataset or a single record). However, in the following, we will use the term *version* as a policy-prescribed reference to a state.

<!--Versioning is intuitive and can trivially be implemented, e.g., by appending version number suffixes to dataset names[^5]. However, since versioning is merely a *policy*, there is no guarantee for enforcement. Further, "There is currently no agreed standard or recommendation among data communities as to why, how, and when data should be versioned" [^6]. The W3C[^7] provides guidance on how revisable data should be handled on the web. However, the (+RDA)[^8] considers these standards too simple and not scalable. 

[^5]: <https://library.stanford.edu/research/data-management-services/data-best-practices/data-versioning>

[^6]: <https://www.ands.org.au/working-with-data/data-management/data-versioning>

[^7]: <https://www.w3.org/TR/dwbp/#dataVersioning>

[^8]: <https://rd-alliance.org/data-versioning-rda-8th-plenary-bof-meeting>
-->

A concern in data versioning is how to reach "[a] consensus about when changes to a dataset should cause it to be considered a different dataset altogether rather than a new version" [^w3org]. This question is futile from a pure identity perspective since every state needs to be identifiable. The distinguishing between version versus state, therefore, is mainly connected to the nature of (+^PID) (and the costs associated with minting them) [@Klump2016], as well as the notion of hierarchical association and provenance.

[^w3org]: [https://www.w3.org/TR/dwbp/#dataVersioning](https://www.w3.org/TR/dwbp/#dataVersioning)

As mentioned above, it is open to debate whether reproducibility is a hard requirement for data citations of a revisable dataset. If not, resolving citations could be deprecated altogether (pessimistic) or only allowed until the data has changed (opportunistic). The opportunistic approach could, e.g., be implemented by timestamping modifications (e.g., by comparing the file systems' last modification date) or through fixity checking. The (+RDA) (+WGDC) [@Rauber2015] elaborates on this approach: A dataset (or a subset) should be given a new identity when the data has changed since the last time the dataset (or subset) was requested. This is recommended to be implemented using a normalized query store and checksums [@Ball2015].

<!-- [@Gray2002] advocate for all non-regeneratable data to remain available forever, and [@Buneman2010] suggests a dataset to remain accessible once it has been cited. A citation can include a version number in their implementation, allowing the system to dereference a citation to the corresponding state. A similar approach is implemented in the dataverse [@Crosas2011], where citations optionally can contain a dataset version number. -->

<!-- A simple way of keeping previous states available is snapshotting, e.g., as implemented by the engineering materials repository MatDB[^10]. However, versioning and snapshotting at fixed points in time do not suit well for (near) real-time data [@Huber2015]. Further, depending on the size and revision frequency, storing all revisions as separate exports may not be feasible [@Rauber2015]. To circumvent these issues [@AltKin07] and the RDA WGDC [@Rauber2015a, Rauber2015, Proll2013] recommend storing state changes on a record level rather than on a dataset level: Every state of every record remains persistently stored and is associated with the CUD operation that touched it.  All CUD operations are timestamped, allowing the identification of every state's validity period. This approach is implemented by [@Alawini2017] for the eagle-iV database.

The concept of versioning at record-level challenges the common understanding of dataset-wide versioning and opens a debate about whether two subsets of two different dataset-level versions containing identical data should share the same identity. 

[^10]: <http://doi.org/10.17616/R3J917>
-->

### Access
There is a common agreement that data citations should use actionable (+^PID) as an access mechanism. E.g., [@AltKin07] suggests a citation to contain an identifier that can be resolved to a landing page (not to the data itself), a requirement also specified by the (+JDDCP) [@JDDCP2014; @Altman2015]. The landing page, in turn, should contain a link to the data resource. The advantage is that the identifier can be resolved regardless of whether the data is behind a paywall or does not exist anymore (at all or in the referenced state).

<!-- <abbr>DOIs</abbr> are a commonly used actionable <abbr>PID</abbr> in data citations, which can be explained by its maturity and the availability of a global resolving service. [@Honor2016]. Alternatives to <abbr>DOI</abbr> are <abbr>ARK</abbr>, <abbr>PURL</abbr>, <abbr>URL</abbr>/permalinks, or <abbr>PURL</abbr> [@Klump2016, Starr2015]. -->

The question of data access is connected to reproducibility and on what granularity level changing states of a revisable dataset should be stored.

### Citation texts
Data citation systems should use metadata standards [@CODATA2013] and be capable of generating human-readable citation strings to facilitate the use of data citations and lower the boundaries for data citation [@Buneman2016; @Rauber2015]. An advantage of citation texts, including, e.g., title, author, and date, is that they allow the reader to quickly assess the relevance, quality, and concurrency of the cited material [@Buneman2010]. The ability to automatically create citation texts is also recommendation 11 for data citations of the (+RDA) (+WGDC).

<!--
The requirements for the metadata attributes to construct human-readable strings slightly vary between recommendations and implementations. An overview of some of the differences can be found in table <a href=" #tab_snippet" data-reference-type= "ref" data-reference= "tab_snippet">1</a>. All the listed recommendations[^11] require the inclusion of Author, Title, and Version. Asides from the stated, the following attributes may also be required: Publish Date, distributor, subset definition, Institution, Related Links, File type, Type of Data, Licences, Project name, Keywords, Repository Name, location (if no resolvable identity is used).

|             | PID  | Subject | Description | UNF  |
| :---------- | :--- | :------ | :---------- | :--- |
| DataCite    | DOI  | Yes     | Yes         | No   |
| Dublin Core | Yes  | Yes     | Yes         | No   |
| Mendeley    | DOI  | No      | Yes         | No   |
| Figshare    | DOI  | Yes     | Yes         | No   |
| Dryad       | DOI  | Yes     | No          | No   |
| Dataverse   | DOI  | No      | No          | Yes  |
| Zenodo      | DOI  | No      | No          | No   |
| RDA WGDC    | Yes  | No      | No          | No   |


[^11]: DataCite: <https://schema.datacite.org>; Dublin Core <https://schema.datacite.org>; Mendeley <https://data.mendeley.com>; Figshare <https://figshare.com/>; Zenodo <https://zenodo.org>; WGDC [@Rauber2015]
-->


# Implementation
The (+OCCUR) enables automated citation creation and dereferencing for data served by (+OPeNDAP) servers. It is implemented as a web service that allows users to:

1. Assign and store identities to data 
1. Create identifiers for identities 
1. Resolve identifiers while verifying that the data has not changed since the creation of the identifier. 
1. Broker identities, i.e., ensuring that identical data shares the same identity.
1. Generate formatted citation snippets for (+OCCUR) identifiers and any arbitrary (+OPeNDAP) query. 

(+OCCUR) is built as a third-party web service without needing to modify data repositories.


## OPeNDAP Primer
The (+OPeNDAP)[^opendap] simplifies access to remote data. It is widely used (but not limited to) to access remote (+HDF) and (+NetCDF) files. An (+OPeNDAP) server allows an (+OPeNDAP) client to request the structure and metadata of data and query subsets of the data. 
The client retrieves those requests through (+HTTP) `GET` requests to an (+OPeNDAP) (+URL) served by the (+OPeNDAP) server. The structure and metadata are retrieved by requesting the (+DDS) and the (+DAS) of the data. The former describes the shape of the data (such as the dimensions), while the latter contains metadata populated by the data provider. Queries for data subsetting are specified through constraint expressions in the (+URL) query strings.

[^opendap]: here, referring to the protocol, not to the company of the same name developing the protocol


## Retrieving identities
In (+OCCUR), data identity is described by:

a) the query (i.e., the full (+OPeNDAP) (+URL)) that produced the data, and 
b) a (+UNF) of the result set of the query (i.e., the data). 

Though not strictly necessary, the time of the identity creation (which serves as a proxy for the query execution time) is used as a third attribute to describe an identity to increase convenience and human readability.

An identity's identifier is the concatenation of the (+OPeNDAP) (+URL) and the identity creation timestamp. (See class diagram in Figure \ref{fig_store}). Identities and their corresponding identifiers are permanently stored within (+OCCUR) upon user request.

A user can retrieve a data identifier by submitting the (+OPeNDAP) (+URL) used to produce the data to (+OCCUR) with the following (+REST) (+API) call:

```http
GET $OCCUR_HOST/store/?dap_url=$DAP_URL
```

When a user requests to retrieve the identity of data, (+OCCUR) first fetches the queried data from the (+OPeNDAP) server and creates a (+UNF) (by default, an MD5[^mdhash] hash) of the current state of the data. (+OCCUR) will then verify if an identity to the same query has already been stored. If not, it stores the query together with the (+UNF) and the current timestamp. The user is then provided with the newly created identifier. 

[^mdhash]: MD5 is a cryptographically broken algorithm. However, we here utilize it merely to verify data integrity (i.e., to identify unintentional corruption). The hash function may be replaced in the future, e.g., by SHA-256, to verify against malicious data corruption. 

If one or more identifiers to the same (+OPeNDAP) (+URL) already have been stored, (+OCCUR) compares the (+UNF) of the current data state to the (+^UNF) of the stored identities. If an identity with the identical query and (+UNF) is found, it is assumed that the data is currently in the same state as when this stored identity was created. Their identity is, therefore, identical, and the user is provided with the identifier of this stored identity. 

If the (+UNF) of the current state differs from all stored (+^UNF) to the same (+OPeNDAP) (+URL), a new identity is created, and the user is provided with the newly created identifier. Figure \ref{fig_store} schematically illustrates the flow for retrieving an identity.

![Class diagram of the data identity and flow for retrieving an identity. \label{fig_store}](images/C1/store.png){ width=90% }


## Dereferencing identities
(+OCCUR) follows an opportunistic approach for dereferencing identities: An identifier can only be dereferenced during the time interval between identity creation and the time the state of the referenced data changes. A user can dereference an identifier by submitting the following (+REST) (+API) call:

```http
GET $HOST/dereference/?identifier=$IDENTIFIER
```

(+OCCUR) first verifies if an identity for the submitted identifier exists/is stored. If so, (+OCCUR) will look up the associated (+OPeNDAP) query and (+UNF). (+OCCUR) then executes the stored query to retrieve the current state of the data. For this state, the (+UNF) is calculated and compared to the stored (+UNF). If the (+^UNF) match, it can be assumed that the data is in the same state as it was during the identity creation. The user thus can be redirected to the (landing page of the) data.

If the (+^UNF) don't match, the current state can be assumed to differ from the state of the data at the identity creation. The user, thus, will be redirected to a landing page containing a warning indicating that the referenced state of the data is not accessible anymore. It then is up to the user whether or not to retrieve the data in its current state. The flow of the dereferencing is illustrated in figure \ref{fig_dereference}.

![Flow for the de-referencing an identity. \label{fig_dereference}](images/C1/dereference.png)


## Formatting citations

(+OCCUR) allows the creation of human-readable formatted citation snippets from both (+OCCUR) identifiers and plain (+OPeNDAP) (+^URL). This service can be understood as the (+OPeNDAP) pendant to [www.crosscite.org](www.crosscite.org), which allows the creation of citation snippets from (+^DOI).

A user specifies the (e.g., journal) formatting style[^csl] together with the identifier or the (+OPeNDAP) (+URL). Besides human-readable strings, a user can choose to retrieve the raw bibliographic information in BibTeX or (+CSL) (+JSON) format. In its simplest form, a user can request a formatted citation snippet with the following (+REST) (+API) calls:

[^csl]: (+OCCUR) supports any of the formats defined in the official repository for (-CSL) citation styles [https://github.com/citation-style-language/styles](https://github.com/citation-style-language/styles).

To format an identifier:

```http
GET $HOST/format/identifier=$IDENTIFIER&style=$STYLE
```

Or, for a plain (+OPeNDAP) (+URL):

```http
GET $HOST/format/dap_url=$DAP_URL&style=$STYLE
```

For the first case, (+OCCUR) will first verify if an identity with the specified identifier exists. If it exists, (+OCCUR) will look up the according (+OPeNDAP) query and derive the (+URL) of the according (+DAS). The (+DAS) is then fetched, and (+OCCUR) will extract metadata from its "global" section[^das_global]. (+OCCUR) will hereby extract any legal citeproc-py keyword (see appendix \ref{app_citeproc}). This metadata is converted into the citeproc-csl (+JSON) format. In case the (+DAS) contains a (+DOI), (+OCCUR) will query <https://www.doi.org> for the (+CSL) (+JSON) representation of this (+DOI)[^content_negotiation] to retrieve additional metadata.

[^das_global]: the part that satisfies: ```/(? <= global.* (?=)/*ims*```

[^content_negotiation]: (+OCCUR) uses content negotiation to query doi.org with the header ```Accept: application/vnd.citationstyles.csl+json``` to retrieve a CSL-JSON metadata response

The (+CSL) (+JSON) of the (+DAS) and the (+CSL) (+JSON) of the (+DOI) are then merged. Together with the style definition, they are fed into citeproc[^citeproc] to produce a formatted citation string, which is then served to the user.

[^citeproc]: (+OCCUR) uses citeproc-py (<https://github.com/brechtm/citeproc-py>)

In case the (+DAS) does not contain a (+DOI) or the user wants to overwrite the (+DOI) specified in the (+DAS), the user can choose to add a (+DOI) to the snippet request, i.e., with the following (+REST) (+API) call:

```http
GET $HOST/format/identifier=$IDENTIFIER&style=$STYLE&doi=$DOI
```

The flow for creating the citation snippet is illustrated in figure \ref{fig_format}.

![Formatting of citations. \label{fig_format}](images/C1/format.png)

## Web frontend and landing pages
The web-service endpoints default to resolve to landing pages rather than to the identifiers or data itself. The reference/identifier endpoint resolves to a landing page that allows the user to inspect the stored identity (specified by the (+URL), the (+UNF), and identity creation timestamp), access the referenced data, as well as to format the identity to a human-readable citation snippet. 

Seminally, the dereference/data endpoint resolves to a landing page of the data rather than the data itself. Those landing pages provide the user with links to the identifier landing page and the referenced data. If the referenced state of the data is no longer available, the landing page additionally contains a warning informing the user that the data has changed. The user may now choose to access the newer state of the data or create a new identity.


## Use Example
A schematic timeline for using (+OCCUR) might look as follows: User A queried data from an  (+OPeNDAP) server. The user now wants to create an identifier for this data to include in a reference. They thus take the (+OPeNDAP) (+URL) used to query the data to (+OCCUR) and request an identifier. (+OCCUR) resolves the (+OPeNDAP) URL, fetches the data, computes the (+UNF), and stores it together with the (+OPeNDAP) (+URL) and the current timestamp as a new identity. The user then is provided with the identities' identifier. They then may use this identifier to create a human-readable citation string. Later, another user, B, may receive the reference, including the (+OCCUR) identifier, and uses (+OCCUR) to dereference it to the referenced data. (+OCCUR) will first verify if the referenced state is still available. If so, the user is provided with the data. If not, the user is given a warning. They now have the option to either access the newer state of the data or create a new identity.

Both the retrieval of an identity and the dereferencing of the identifier are illustrated in figure \ref{fig_storingcitation}.

![Timeline of identity creation and retrieval. \label{fig_storingcitation}](images/C1/storing_resolving.jpg)

# Related work
<!-- Not obvious if this section is needed considering the detailed description in Rauber2021 -->

The following is a non-comprehensive summary of services for automized data citation creation that have been presented in the past:

MatDB[^MatDB] implements a data publication and citation service for engineering materials. Datasets are made citable by enforcing minimal discipline-specific metadata and minting DataCite (+^DOI). Fixity is assured by snapshotting the dataset at the time of (+DOI) minting. Revisability of data is made possible through policy-enforced versioning [@Austin2016].

[^MatDB]: [http://doi.org/10.17616/R3J917](http://doi.org/10.17616/R3J917)

[@Alawini2017] created a citation service for the (+RDF) eagle-i database. Since this database itself does not version its data (only the most recent version is available), the authors implemented an external service that versions eagle-i data to provide users access to revised data. The service tracks and stores every change in the original dataset. The authors note that this approach is viable for eagle-i since the dataset changes very slowly.

The dataverse networks software [@Crosas2011] aggregates data in "studies." Studies may contain several datasets, and each study shares a common persistent identifier. A citation to a dataset (and subsets) is implemented as the combination of the studies' (+PID) appended with the (+UNF) of the cited data.

[@Cook2016] presents the data product citations at a (+ORNL) (+DAAC). The (+DAAC) assigns a (+DOI) per dataset, which may contain between one and tens of thousands of files. A single file within a dataset can be identified by appending the file's (+UUID) to the (+DOI) (using the ```urlappend``` functionality of the (+DOI) resolver). The (+DAAC) also provides a (+MODIS) subsetting service. The user can request a citation to the subset, comprised of the dataset's citation appended with a textual description of the temporal and/or spatial subsetting. The citation can be requested as a formatted string and in BibTex format.

A very similar approach is implemented for the <abbr>ARM</abbr> Data Archive [@Prakash2016]. Upon data order fulfillment, the user is provided with both the data and a citation that contains a citation, including a textual description of the temporal and/or spatial subsetting. The (+ARM) Data Archive also hosts a citation creator, a (+GUI), that allows the creation of a data citation subject to a fully qualified dataset stream name and optionally manually specified subsetting parameters. The user hereby can choose between the custom (+ARM) citation style, (+APA), (+MLA), and Chicago.

[@Honor2016] describes a reference implementation for data citations of a database holding neuroimaging (the specific use case is the (+NITRC)). The citable units for their use-case are individual images aggregated in studies/projects. Upon data upload, hierarchically, each image and each study is assigned a (+DOI). The authors recognize that this implementation may result in the generation of many (+^DOI); however, they evaluate the solution feasible for their use case.

[@Proll2013] present a reference implementation for data citations to data managed by a (+RDBMS). The system is based on the premise that a timestamped SELECT query can correctly identify data. To allow revisability, the system timestamps every (+CUD) operation and acknowledges the validity time ranges of records rather than allowing modification of records. When a user wants to create a citable subset, the system will store the according timestamped `SELECT` query, calculate a hash for the result set, and assign a (+PID) to the query.

[@Buneman2010] describes an approach for citing digital archives described by the (+EAD) and for the (+IUPHAR) database. They hereby focus on the concept of citable units and their application on hierarchical/structured and revisable datasets.

[@Rauber2021] details the experiences of several adapters of the 14 (+RDA) (+WGDC) recommendations. Among those are the (+CBMI), (+VAMDC), (+CCCA), and the (+EODC). All adapters were able to improve their data distribution infrastructure through the adaptations of the recommendations. However, each adapter's individual implementations and the required work vastly varied. A common challenge can be identified as enabling data versioning and timestamping.

Apart from those fully integrated data citation systems, we want to mention the crosscite (+DOI) Citation formatter[^crosscite]. The crosscite (+DOI) Citation formatter generates citation strings from metadata retrieved from (+DOI) landing pages. The citation strings are formatted through citeproc subject to styles defined through the (+CSL)[^citationstyles].

[^crosscite]: [https://citation.crosscite.org/](https://citation.crosscite.org/)

[^citationstyles]: <https://citationstyles.org/>

# Discussion and Outlook

We demonstrated how a third-party data citation system can be built for data accessible through the RESTful data access protocol (+OPeNDAP). We define a _citeable unit_ as any result set of a `SELECT`/`GET` request. Identities and the corresponding identifiers to citable units can be created and retrieved upon user request. The identities are defined and permanently stored as a combination of the query, the result set's hash, and the identity creation timestamp. Identifiers can be opportunistically dereferenced as long as the referenced identity is still available. The availability is verified by comparing the hash of the referenced identity and the hash of the query result set at the time of dereferencing. It is hereby irrelevant if the data changed intentionally subject to a revision or unintentionally due to rot or malicious manipulation. It shall be noted that the used hashes stay opaque to the user as (+OCCUR) handles the result set verification internally.

(+OPeNDAP) allows users to request data in a delivery container format that differs from the storage container format. (+OCCUR) calculates hashes in the delivery container rather than in the storage format. Regarding hashes, (+OCCUR) requires fetching the complete result set to calculate the hashes, which may inflict high transfers on both (+OCCUR) and the (+OPeNDAP) resource. To avoid these high transfers, hashes of result sets should be calculated at the point of storage and exposed/queried through an (+API). DAP4 already provides the option of including CRC32 checksums in the (+DMR)[^dmr], which could be used for this purpose.

[^dmr]: [https://docs.opendap.org/index.php/DAP4:_Specification_Volume_1#Checksums](https://docs.opendap.org/index.php/DAP4:_Specification_Volume_1#Checksums)

Identities and identifiers are created, stored, and resolved within OCCUR. This comes with the advantage of low cost. However, these identifiers cannot be considered persistent. It is necessary to use an external service such as (+DOI), identifiers.org[^identifiersorg], or name2thing[^name2thing] to mint persistent identifiers to (+OCCUR) identities.

[^identifiersorg]: [https://identifiers.org/](https://identifiers.org/)

[^name2thing]: [https://n2t.net/](https://n2t.net/)

(+OCCUR) identifiers can be converted into human-readable citation snippets as well as machine-readable bibliographic entries in BibTex, (+RIS), and (+CSL) - (+JSON) format. (+OCCUR) hereby exploits that (+OPeNDAP) resources have a designated location for dataset-level metadata: the (+DAS). The metadata extracted from the DAS is combined with the hash, the (+OPeNDAP) query literal, and the identity creation time to create snippets through citeproc. We recognize that the (+DAS) might be a duplicate metadata location since a dataset-level (+DOI) may have already been registered in many cases. In these cases, we encourage merely the inclusion of the dataset (+DOI) into the (+DAS). (+OCCUR) will parse this (+DOI) and resolve it for the corresponding metadata.

(+OCCUR) is meant to be a demonstration that can be adapted to other RESTful data access services such as (+WMS) and (+WCS). It could be envisaged to remain a centralized third-party system or be integrated into the data services.

# Acknowledgments {.unlisted .unnumbered #c1_acknowledgments}
The development of (+OCCUR) was supported by an (+NSF) grant (Award Number: 1302212)[^nsfgrant] and the (+ESIP) federation as a 2018 (+ESIP) lab project. We would like to thank for the feedback of we received at the 2017 Bren Ph.D. Symposium and the 2018 (+ESIP) summer meeting; particularly James Gallagher, Ted Habermann, and Mark Parsons.

[^nsfgrant]: [https://www.nsf.gov/awardsearch/showAward?AWD_ID=1302212](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1302212)

# Competing interests {.unlisted .unnumbered}

The authors have no competing interests to declare.

\clearpage

# Appendix {-}

## RDA WGDC recommendation compliance
OCCUR implements several of the 14 (+RDA) (+WGDC) recommendations [@Rauber2015]. It stores queries with metadata, (+PID), and timestamps (R7, R8, R9). It also uses (+^UNF) to achieve result set verification (R6). Since (+OPeNDAP) normalizes queries by alphabetically sorting the constraint expressions, (+OCCUR) can also fulfill the recommendation to ensure query uniqueness and that a single subset is not referred to by more than one identifier (R4). Further, when accessing files through (+OPeNDAP), a stable sorting of the result sets can be guaranteed (R5). (+OCCUR) is machine actionable through its (+REST) (+API) (R12) and resolves to human-readable landing pages (R11), both of which allow the creation of formatted citation snippets (R10). Since (+OCCUR) is a third-party service, it does not address data versioning (R1) and timestamping (R2), as we see those as the repositories' responsibility. Additionally, a technology migration (R13) and migration verification (R14) will only be possible to be carried out in collaboration with the according repository.

Figure \ref{WGDC} evaluates the compliance of the current state of (+OCCUR) with these recommendations. 

![Evaluation of WGDC compliance. \label{WGDC}](images/C1/comparison.png)

## Legal citeproc keywords
\label{app_citeproc}

```python
NAMES = ['author', 'collection_editor', 'composer', 
		 'container_author', 'editor', 'editorial_director', 
		 'illustrator', 'interviewer', 'original_author', 
		 'recipient', 'translator']
DATES = ['accessed', 'container', 'event_date', 'issued',
         'original_date', 'submitted']
NUMBERS = ['chapter_number', 'collection_number', 
           'edition', 'issue', 'number', 'number_of_pages', 
           'number_of_volumes', 'volume']
VARIABLES = ['abstract', 'annote', 'archive', 
			 'archive_location', 'archive_place', 
			 'authority', 'call_number', 'citation_label', 
			 'citation_number', 'collection_title', 
			 'container_title', 'container_title_short', 
			 'dimensions', 'DOI', 'event', 'event_place',
			 'first_reference_note_number', 'genre',
             'ISBN', 'ISSN', 'jurisdiction', 'keyword', 
             'language', 'locator', 'medium', 'note', 
             'original_publisher', 'original_publisher_place', 
             'original_title', 'page', 'page_first', 'PMCID', 
             'PMID', 'publisher', 'publisher_place', 
             'references', 'section', 'source', 'status', 
             'title', 'title_short', 'URL', 'version', 
             'year_suffix']
```




# References{-}
