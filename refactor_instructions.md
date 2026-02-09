# Task

Change the current implementation of the aligner program to a system where users designate forward primers, reverse primers and probes seperately. The tool should then specifically search for an amplificate between forward and reverse primers. It should match the probe in the amplificate region. It should then generate a similar report as currently (excel export being the most important output), which specifically analyses the minimum number of missmatches among the best fitting primers and probes.

## UI

### Sequences

Currently, users choose a sequence file as fasta. They then either import an oligo file as fasta or copy paste a fasta formatted list in the text field. There is currently no distinction between oligo types. There is also still an "output file" function.

Change this into the following:

1. Keep the "Sequence file" system unchanged. Users will continue to import their reference sequence alignments as fasta.
2. Remove the "oligo file" functionality
3. Remove the "output file" functionality in the "input files" area. This will be a functionality for the report dialogue only.
4. For the "Oligo sequences" area, make the following changes:
  - Implement buttons for "Forward primers", "Reverse primers" and "Probes"
  - Each should be their own bucket of oligo sequences. Clicking one will clear the current content of the box and load the previously saved sequences for that category (if there are any)
  - Buttons should be red if no valid fasta sequences are saved and green if valid fasta sequences are saved.
  - Provide a "save" and "load" button. upon saving, the program should generate a json file that contains all saved oligos, split into the three categories. Loading such a json should fill the relevant categories automatically. 
  - Running an analysis should only be possible if there is at least 1 valid forward primer and 1 valid reverse primer.

### Quick settings area

Keep the same but "min Oligos matched" shoudl now be per category. Ambiguity match display can default to "Dots", but have letters as an option.

### Amplicon constraints

This can remain unchanged but has to harmonize with the new analysis methodology.

### Buttons

Keep "Run analysis". Have this in red while there is not at least one forward and one reverse primer saved. Eliminate the "File only" and "Export excel" buttons in that area.

## Analysis method

Currently, all oligo sequences are treated the same. The program searches for possible amplicons within a certain length range and then performs the aligment for assembling the different patterns.

Change this to the following:

1. Amplificates should only be possible between designated forward and reverse primers. Depending on the orientation of the input, reverse primers can be sense orientation and forward in reverse orientation.
2. Probes have to be matched between at least 1 forward and 1 reverse primer to be valid for the analysis.
3. If primers can be matched within the length constrains, evaluate all possible forward and reverse primers that can be matched. (similar to how it is implemented currently)

*--> Intention behind these changes:
The idea is that for a given qPCR assay, only one well matching forward primer, reverse primer and probe are needed to generate robust signals. The tool has to be able to show if at least one such combination exists for any given sequence variation. Variations where there are no combinations without missmatches or no combinations without substantial numbers of missmatches, these are an inclusivity concern for the assay* 

## Output

The current implementation of the output is already well suited for purpose. Keep the presentation mostly the same but harmonize with the altered analysis method. The excel-file is the most important output. 

For the output, organize the oligos in the following sequence: (left to right) 
1. Forward primers
2. Probes
3. Reverse primers

Clearly indicate where one section ends and where the other starts. Within category, sort them by input sequence
