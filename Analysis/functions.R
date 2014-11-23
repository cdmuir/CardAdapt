exportTable <- function(file, data, colnames, rownames, h1 = "", ...)
{
	# Check that there are the right number of columns and rows
	if(length(colnames) != ncol(data)) stop("number of columns in data does not match length of column names")

	if(length(rownames) != nrow(data)) stop("number of rows in data does not match length of row names")

	table <- file(file, open = "w")
	cat(h1, "&", colnames[1], file = table)
	if (length(colnames) > 1) 
	{
		for (i in 2:length(colnames)) cat(" &", colnames[i], file = table)
	}
	cat(" \\\\ \\midrule\n", file = table)
	
	for (i in 1:length(rownames))
	{
		cat(rownames[i], "&", data[i, 1], file = table)
		if (length(colnames) > 1) 
		{
			for (j in 2:length(colnames)) cat(" &", data[i, j], file = table)
		}
		cat(" \\\\\n", file = table)
	}
	
	close(con = table)	
	
}