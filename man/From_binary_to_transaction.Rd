\name{From_binary_to_transaction}
\alias{From_binary_to_transaction}
\title{Creation of data as transaction from binary data}
\usage{
transactiontoBitmax(Binary_data,variable_names,deli)
}
\arguments{
\item{Binary_data}{A Binary Dataframe}

\item{variable_names}{The colnames of you Binary Data}

\item{deli} {The delimitator you wan't to set for your transaction format}
}

\description{
Returns the corresponding Transaction Data.
}

\examples{

data("Transacfruits")
binaryfruits = From_transaction_to_binary (as.matrix(Transacfruits),",")
Transaction_with_point_separator = From_binary_to_transaction(binaryfruits,names(binaryfruits),".")
## Import transaction data Transacfruits,create the correspondant binary matrix and recreate the transaction format with another delimitator
}
