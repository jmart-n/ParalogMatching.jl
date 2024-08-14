using ParalogMatching

function computePairing(firstAlignment::ParalogMatching.Alignment, secondFasta::AbstractString, outputFile::AbstractString)::nothing
    secondAlignment = ParalogMatching.read_fasta_alignment(secondFasta, header_regex=header_regex)
    preparedFastas = ParalogMatching.prepare_alignments(firstAlignment, secondAlignment, cutoff=0)
    matchingResult = ParalogMatching.run_matching(preparedFastas, batch=1) #default pseudocount=0.5, covariation method
    ParalogMatching.write_fasta_match(preparedFastas, matchingResult, outputFile)
end;

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

global header_regex = r"^(?<id>.+)\|TaxID=(?<species>.+)\|\|"

fixedFasta = ARGS[1]
inputFasta = ARGS[2]
outputFolder = ARGS[3]

if !isdir(outputFolder)
    mkdir(outputFolder)
end
inputFastaName = split(inputFasta, "/")[2]
outputFastaName = outputFolder * "/" * join(insert!(split(inputFastaName, "."), 2, "paired"), "_", ".")

fixedAlignment = ParalogMatching.read_fasta_alignment(fixedFasta, header_regex=header_regex)

computePairing(fixedAlignment, inputFasta, outputFastaName)