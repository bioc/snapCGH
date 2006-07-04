
"readPositionalInfo" <-
  function(input, source, path = NULL){
    if(!is.null(path))
      file <- file.path(path, RGList$targets[1])
    else
      file <- input$targets[1]

    source <- match.arg(source, c("agilent", 
        "bluefuse", "nimblegen"))
   
    switch(source,
           agilent = {
             split <- strsplit(input$genes$SystematicName, split = ":")
             chr <- start <- end <- vector(length = length(split))
             for(i in 1:length(split)){
               if(length(split[[i]]) == 1) {
                 chr[i] <- start[i] <- end[i] <- NA
               }
               else {
                 chr[i] <- substr(split[[i]][1], 4, 5)
                 if(chr[i] == "X")
                   chr[i] = 23
                 else if(chr[i] == "Y")
                   chr[i] = 24
                 
                 temp <- strsplit(split[[i]][2], "-")
                 start[i] <- temp[[1]][1]
                 end[i] <- temp[[1]][2]
               }
             }
           },
           bluefuse = {
             if(!is.null(path)){
               n <- nchar(path)
               if(substr(path, n, n) != "/")
                 { path = paste(path, "/", sep = "") }
               file <- paste(path, input$targets$FileName[1], sep = "")
             }
             else {
               file <- input$targets$FileName[1]
             }
             skip <- readGenericHeader(file, columns = c("ROW"))$NHeaderRecords
             table <- read.table(file, skip = skip, header = T, sep = "\t")
             chr <- start <- end <- vector(length = nrow(table))
             for(i in 1:length(table$CHROMOSOME)){
               if(table$CHROMOSOME[i] == "-1"){
                 chr[i] <- start[i] <- end[i] <- NA 
               }
               else {
                 chr[i] <- table$CHROMOSOME[i]
                 start[i] <- end[i] <- table$POSITION[i]
               }
             }
           },
           nimblegen = {
             chr <- as.vector(unlist(input$genes$SEQ_ID))
             chrnum <- sapply(chr, function(z){strsplit(z, "chr")[[1]][2]})
             indx <- which(chrnum == "X")
             indy <- which(chrnum == "Y")
             chr[indx] = 23
             chr[indy] = 24
             input$genes$Chr <- as.numeric(chr)
             names(input)[5] <- "Position"
             input$Position <- as.numeric(input$Position)/1000000
           }
             )
    input$genes <- data.frame(input$genes, Chr = as.numeric(chr), Start = (as.numeric(start)/1000000), End = (as.numeric(end)/1000000))
    input
  }
         
