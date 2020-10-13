% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/average_marker_intensity_within_radius.R
\name{average_marker_intensity_within_radius}
\alias{average_marker_intensity_within_radius}
\title{average_marker_intensity_within_radius}
\usage{
average_marker_intensity_within_radius(
  sce_object,
  reference_marker,
  target_marker,
  radius = 20
)
}
\arguments{
\item{sce_object}{A SingleCellExperiment object in the form of the output of format_image_to_sce}

\item{reference_marker}{A string specifying the marker that is used for reference cells}

\item{target_marker}{A string specifying the marker to calculate its average intensity}

\item{radius}{An integer specifying the radius of search for cells around the reference cells}
}
\value{
A single number is returned
}
\description{
Calculates the average intensity of the target_marker within a radius
from the cells positive for the reference marker.
Note that it pools all cells with the target marker that are within the specific radius of
any reference cell. Results represent the average intensities within a radius,
but do not correspond to metrics for each cell
}
\examples{
average_marker_intensity_within_radius(SPIAT::formatted_image,
                                       reference_marker ="CD8",
                                       target_marker = "CD4",
                                       radius=30)
}