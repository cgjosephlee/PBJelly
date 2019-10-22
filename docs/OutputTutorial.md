# PBJelly Output Tutorial
https://sourceforge.net/p/pb-jelly/wiki/PBJelly%20Output%20Tutorial/

*Below is the response to user Nikki's questions about PBJelly's output. My response text is in **bold**. This is the beginning of a fuller PBJelly Tutorial.*

I'm using PBjelly to improve my genome assembly using long sequence reads.
The run completes successfully, and I'm now trying to understand what has
been changed in the assembly and why. This has resulted in a couple of
questions on the output generated:

1) Why and how is the assignment of the new reference names done by PBJelly?

For example, I have 3 reference scaffolds, that after renaming by PBJelly have names that look like:
```
9250459_1162_24830|ref0158577|ref0122806|ref0096825
4075816_1513_39396|ref0142504|ref0052326|ref0158577
11012802_95505_2694272_561104-,...,10967441+|ref0226469|ref0158577|ref0182998
```
Why are 3 new names assigned to 1 scaffold? (e.g. ref0158577, ref0122806 and ref0096825 to the first scaffold?)

Why are those "subnames" (e.g. ref0158577) not unique?

**During PBJelly’s setup stage, we add a key to the end of every scaffold name with the structure: `ref\d{7}`.\
This is a unique identifier that PBJelly uses in order to avoid special characters which are allowed in Fasta Entry Headers but not in Unix FileNames (or are at least frowned upon) from affecting the temporary files PBJelly creates. When you run multiple iterations of PBJelly, a new `ref\d{7}` is added to the end of the scaffold name. However, these is an option in the Setup stage for “—noRename” that will allow you to use the existing `ref\d{7}` key in the fast entry headers.**

2) In the gap_fill_status.txt file I see that a gap gets filled in ref0158577:
```
ref0158577.1e3_ref0158577.2e5 filled
ref0158577.2e3_ref0158577.3e5 filled
```
What do the .1 and .2 mean in ref0158577.1e3_ref0158577.2e5?

How do I figure out in which of the three original scaffolds (that all have a new name containing 'ref0158577') have been filled?

**For your first question (.1/.2), let’s start from the beginning.\
Let’s say you have a single scaffold with 5 gaps. Setup will name the scaffold ref0000001 and identify the 5 gaps. If there are 5 gaps in your scaffold, that means there are 6 contigs. Gaps are defined by the flanking sequence around them. PBJelly uses python which uses 0-based indexing, So your first contig in your scaffold is 0, second is 1, etc. So your first gap is named ref0000001_0_1. If you look in gapInfo.bed you’ll see these names with coordinates of where your gap is located within the scaffold.\
This is the system I started within, however, once PBJelly was expanded to include inter-scaffold gap-filling, I had to expand on this functionality a little bit. Using the same ideas I just described, we now needed to build a graph to represent the genome AND directionality for inter-scaffold gaps. Nodes in the graph represent contig ends and edges represent sequence. For example, The same ref0000001_0_1 gap would be located between the 3’ end of the 0th contig and the 5’ end of the 1st contig. Similarily, the 1st contig is represented as ref0000001.1e5_ref0000001.1.e3.\
In order to tie this information back to your original scaffold, you’d need to build a hash using gapInfo.bed with key holding the `ref\d{7}` key and value the original name.**

3) In jelly.out.fasta the scaffolds have new names looking like Contig0
etc. Is there somewhere where I can find back from which of the original
scaffolds this contig is derived?

**Check check the lift-over table.json as a start. Soon I’ll be refactoring the lift-over table so it’s more useful. Briefly, you’ll have key:value of new Contig name followed by a list of all the graph edges that comprise the new contig.**

4) In gap_fill_status.txt I see the following types of changes:
```
doubleextend
filled
nofillmetrics
overfilled
singleextend
```
- What does 'overfilled' mean? Is this an inter-scaffold connection? If not, how can I see which inter-scaffold connections have been made (if any)?

- In jelly.out.fasta the smallest contig is \~300 nt, while in my input file the smallest scaffold was 1000 nt. This indicates that pbjelly has split something up. Where can I find back why one or more scaffolds have been split up, and which ones?

- Could you give a short explanation on what a line like this in the gapInfo.bed file means:
```
10515137_3627_104321_1776636+,...,7058545+|ref0197248|ref0108226|ref0155877
na na ref0155877_0_0 3
```

**doubleextend means we extended both flanking sequences into the gap\
singleextend means we extend one of the flanking sequences into the gap\
filled means we were able to connect the flanking sequences and remove the gap\
nofillmetrics means we were unable to create any kind of consensus sequence to fill the gap\
overfilled means we extended one or both flanking sequences into the gap, however the amount of sequence we placed into the gap is greater than the predicted size of the gap (based on the number of Ns) and yet we were unable to find a sequence that spans the gap. For example, if we have a 1000bp gap, and we extend the e3 flanking sequence 800bp and the e5 flanking sequence 1000bp, but the two extending sequences aren’t united, we’ve overfilled the gap by 800bp.**

**PBJelly will trim contigs if it finds an appropriate situation to do so. However, this functionality is slightly finicky and I’m actively developing it and PBHoney in order to address this problem. Briefly, many gaps are the result of one allele being constructed on one flanking sequence, and the other allele being constructed on the other flanking sequence. PBJelly can occasionally identify these regions based on mapping unto the end of one contig and slightly ‘inland’ on the other contig. PBJelly will occasionally trim one of the alleles and try to close the gap. This is something I’m hoping to improve upon, so anything you learn about this kind of behavior will be very beneficial to share.**

**That’s a singleton with no gaps.\
cur start end id start/end_flag\
The start/end_flag can be interpreted as : 0x1 == first contig in scaffold; 0x2 == last contig in scaffold.**

**This is another aspect of how I had to overhaul the previous PBJelly intra-scaffold only gap-filling to inter-scaffold gap filling. A line such as the above helps when building the graph that represents the genome. gapInfo.bed represents all of the gaps in the genome, including inter-scaffold gaps. For singleton contigs, it has two gaps it needs to represent - the gap going from it’s 5’ and 3’ ends.**

Sorry for all the questions, and thanks a lot for your help!

With kind regards,\
Nikkie
