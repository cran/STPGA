\name{makeonecross}
\alias{makeonecross}

\title{Make a cross from two solutions and mutate. }
\description{
Given two training sets, identifiers for candidates, this function makes a new solution using crossover and one point mutation with probability mutprob. Only one mutation is allowed. 
}
\usage{
makeonecross(x1, x2, Candidates, mutprob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
a vector of identifiers selected from the candidate set.
}
  \item{x2}{
a vector of identifiers selected from the candidate set that has the same length as x1.
}
  \item{Candidates}{
vector of identifiers for individuals in the candidate set. 
}
  \item{mutprob}{
point mutation probability for each individual generated. Only one mutation per solution is allowed. 
}
}

\author{
Deniz Akdemir
}
