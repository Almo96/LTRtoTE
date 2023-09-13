import os

def controlIO(LTR, genome, output, output_path, minlenTE, maxlenTE, tsd):
    # Check if the extensions are .fasta or .fa
    valid_extensions = ['.fasta', '.fa']
    for arg_name, arg_value in [("LTR", LTR), ("genome", genome)]:
        if not arg_value.endswith(tuple(valid_extensions)):
            raise ValueError(f"{arg_name} must have a .fasta or .fa extension.")


    print("LTR:", os.path.basename(LTR))
    print("Genome:", os.path.basename(genome))
    print("Output path:", output_path)
    print("Output:", output)
    print("Minimum length:", minlenTE)
    print("Maximum length:", maxlenTE)
    print("TSD:", tsd)
