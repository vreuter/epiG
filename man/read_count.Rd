% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiG-fetch.R
\name{read_count}
\alias{read_count}
\title{Count Reads}
\usage{
read_count(file, refname, start, end)
}
\arguments{
\item{file}{path to bam file}

\item{refname}{reference name}

\item{start}{start of region}

\item{end}{end of region}
}
\value{
number of reads and total bps in reads overlapping region
}
\description{
Count reads overlapping the specified region
}
\examples{
# Retrieve paths to raw data files
bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")

read_count(bam_file, "chr20", 57400000, 57400000 + 100)
}
\author{
Martin Vincent
}

