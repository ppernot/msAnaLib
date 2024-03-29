% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compressMS.R
\name{compressMS}
\alias{compressMS}
\title{Down-sample a FT-ICR MS by choosing points closest to a preset
regular grid or around a set of targets.}
\usage{
compressMS(
  file_in,
  file_out,
  ms_type = "fticr",
  mzMin = 70,
  mzMax = 250,
  dmz = 0.002,
  test = FALSE,
  compMode = "grid+targets",
  targets = NA,
  dmzTarget = 0.5
)
}
\arguments{
\item{file_in}{(string; mandatory) Input DMS file to be compressed}

\item{file_out}{(string; mandatory) Output (compressed) file}

\item{ms_type}{(string; optional) Type of mass spectrum,
one of `esquire` or `fticr` (default).
Compression only for `fticr` MS.}

\item{mzMin}{(numeric; optional) Min m/z value of compression grid}

\item{mzMax}{(numeric; optional) Max m/z value of compression grid}

\item{dmz}{(numeric; optional) m/z step or compression grid}

\item{test}{(logical; optional) if TRUE, shodt run to check the
arguments validity and the compression efficiency}

\item{compMode}{(string; optional) Compression mode:
either `grid`, `targets` (default),
or their combination `grid+targets`.}

\item{targets}{(string; optional)  Path to targets file.}

\item{dmzTarget}{(numeric; optional) Half-width of m/z interval
centered on target reference m/z value.}
}
\value{
Returns nothing. It is used for its side effects.
}
\description{
Down-sample a FT-ICR MS by choosing points closest to a preset
regular grid or around a set of targets.
}
