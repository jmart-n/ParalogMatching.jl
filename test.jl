function parseInputs(inputFolder::AbstractString, outputFolder::AbstractString)::Tuple{Vector,Vector}
    # Check/create directories
    if !isdir(inputFolder)
        println("Check input folder exists")
        exit(0)
    end
    if !isdir(outputFolder)
        mkdir(outputFolder)
    end

    # build lists
    inputFiles = readdir(inputFolder)
    outputFilenames = [outputFolder * "/" * join(insert!(split(filename, "."), 2, "paired"), "_", ".") for filename in inputFiles]
    inputFilenames = [inputFolder * "/" * filename for filename in inputFiles]
    return inputFilenames, outputFilenames
end;

inputfiles, outputfiles = parseInputs("test", "doubletest")
for (inputFasta, outputFasta) in zip(inputfiles, outputfiles)
    println(inputFasta, outputFasta)
end

