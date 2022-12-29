from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from rich.table import Table
from rich.console import Console
from rich import box
import re
from rich.prompt import Prompt

#User gives the sequence they want to amplify. Can be any length, must include only A, G, C, and T for GC clamp/melting temp, etc. to work properly.
#FASTA = input('\r\nWelcome to Basic Primer Designer. \r\n\r\n** The sequence must only contain A, G, C, and T nucleotides for certain features to work. \r\n** The sequence must be one continous line. Ensure there are no spaces, symbols, or indents (eg. ¶). \r\n\r\nPlease paste the complete 5\'-3\' FASTA sequence you want to amplify, then press enter:')
FASTA = Prompt.ask("[yellow]Welcome to Basic Primer Designer![/yellow] \r\n\r\n[blue]** The sequence must only contain A, G, C, and T nucleotides for certain features to work. \r\n** The sequence must be one continous line. Ensure there are no spaces, symbols, or indents (eg. ¶).[/blue] \r\n\r\n[yellow]Please paste the complete 5\'-3\' FASTA sequence you want to amplify, then press enter:[/yellow]")
my_seq = Seq(FASTA)
str_seq = str(my_seq)

#User specifies basic primer length. For the GC clamp primer, this sets the min bp length but not the max.
primer_sizei = Prompt.ask('\r\n[yellow]What length (in bp) do you want your primer to be (not including res. sites)(press enter for the reccomended 10 bp primer size):[/yellow]', default=10)
try:
    primer_size = round(int(primer_sizei))
except: 
    print('[red]Error: Invalid primer length given. Please try again with an integer response.[/red]')


#nucleotides for res. enzyme site.
res_site_extra = 'AAAAAA'

#predefined restriction enzymes for ease-of-use: forward primer.
res_sitefw = Prompt.ask('\r\nWhat restriction site do you want for your FORWARD primer? (If no site, press enter):')
if len(res_sitefw) < 1:
    res_sitefw = ''
    res_site_extra = ''
elif res_sitefw.lower() == 'hind111':
    res_sitefw = 'AAGCTT'
elif res_sitefw.lower() == 'xho1':
    res_sitefw = 'CTCGAG'
elif res_sitefw.lower() == 'nhe1':
    res_sitefw = 'GCTAGC'
elif res_sitefw.lower() == 'sal1':
    res_sitefw = 'GTCGAC'
else:
    res_sitefw = Prompt.ask('\r\n[red]***Error, this restriction enzyme is not in the database yet. [/red][yellow]Please type the 5\'-3\' recognition site, then press enter (no symbols or spaces):[/yellow]')

#predefined restriction enzymes for ease-of-use: reverse primer.
res_siterv = input('\r\nWhat restriction site do you want for your REVERSE primer? (If no site, press enter):')
if len(res_siterv) < 1:
    res_siterv = ''
    res_site_extra = ''
elif res_siterv.lower() == 'hind111':
    res_siterv = 'AAGCTT'
elif res_siterv.lower() == 'xho1':
    res_siterv = 'CTCGAG'
elif res_siterv.lower() == 'nhe1':
    res_siterv = 'GCTAGC'
elif res_siterv.lower() == 'sal1':
    res_siterv = 'GTCGAC'
else:
    res_siterv = input ('\r\n***Error, this restriction enzyme is not in the database yet. Please type the 5\'-3\' recognition site, then press enter (no symbols or spaces):')

#GC clamp forward primer
#30 nucleotides max GC primer length
gc_fw_1 = my_seq[0:primer_size]
gc_fw_2 = gc_fw_1.find("C", primer_size, primer_size + 10)
gc_fw_3 = gc_fw_1.find("G", primer_size, primer_size + 10)

 #Finds the first G and C within the specified bounds, uses G or C depending on which comes first.
 # If none are found, GC primer will be blank.
if gc_fw_2 < 0:
    gcf_end = -1
    res_sitefw = ''
    res_site_extra = ''
elif gc_fw_3 < 0:
    gcf_end = -1
    res_sitefw = ''
    res_site_extra = ''

elif gc_fw_2 < gc_fw_3:
    gcf_end = gc_fw_2
elif gc_fw_3 < gc_fw_2:
    gcf_end = gc_fw_3
    
else:
    gcf_end = -1

gc_fw_end = gc_fw_1[0:gcf_end + 1]

#Reverses FASTA, gets the first (technically the last in the original FASTA) 30 nucleotides, reverses again, reverse + complements.
#This circumvents python from splicing the last nucleotide index when using [0:-1].
#Biopython lacks a .reverse() so .complement() has to be used twice.
gcrv_1 = my_seq.reverse_complement()
gcrv_2 = gcrv_1.complement()
gcrv_seq1 = gcrv_2[0:30]
gcrv_3 = gcrv_seq1.reverse_complement()
gcrv_4 = gcrv_3.complement()
gcrv_5 = gcrv_4.reverse_complement()

gcrv_seq = gcrv_5

#GC clamp reverse primer
#30 nucleotides max GC primer length
gc_rv_1 = gcrv_seq
gc_rv_2 = gc_rv_1.find("C", primer_size, primer_size + 10)
gc_rv_3 = gc_rv_1.find("G", primer_size, primer_size + 10)

#Finds the first G and C within the specified bounds, uses G or C depending on which comes first.
# If none are found, GC primer will be blank.
if gc_rv_2 < 0:
    gcr_end = -1
    res_siterv = ''
    res_site_extra = ''
elif gc_rv_3 < 0:
    gcr_end = -1
    res_siterv = ''
    res_site_extra = ''
elif gc_rv_2 < gc_rv_3:
    gcr_end = gc_rv_2
  
elif gc_rv_3 < gc_rv_2:
    gcr_end = gc_rv_3

else:
    gcr_end = -1

gc_rv_end = gc_rv_1[0:gcr_end + 1]

#Basic forward primer
fw_seq = my_seq[0:primer_size]

#Basic reverse primer
#Reverses FASTA, gets the first (technically the last in the original FASTA) 30 nucleotides, reverses again, reverse + complements.
#This circumvents python from splicing the last nucleotide index when using [0:-1].
#Biopython lacks a .reverse() so .complement() has to be used twice.
rv_1 = my_seq.reverse_complement()
rv_2 = rv_1.complement()
rv_seq1 = rv_2[0:primer_size]
rv_3 = rv_seq1.reverse_complement()
rv_4 = rv_3.complement()
rv_5 = rv_4.reverse_complement()

rv_seq = rv_5

#Defines the "total" primer including any res. sites.
forward_primer = res_site_extra + res_sitefw + fw_seq
reverse_primer = res_site_extra + res_siterv + rv_seq

gc_forward_primer = res_site_extra + res_sitefw + gc_fw_end
gc_reverse_primer = res_site_extra + res_siterv + gc_rv_end

#Defines the sequence that is amplified including any res. sites.
amplicon = res_site_extra + res_sitefw + str_seq + res_site_extra + res_siterv

#Seperates all A, G, C, and T nucleotides from each primer
mfp_c = re.findall("C", str(forward_primer))
mfp_g = re.findall("G", str(forward_primer))
mfp_a = re.findall("A", str(forward_primer))
mfp_t = re.findall("T", str(forward_primer))

mgcfp_c = re.findall("C", str(gc_forward_primer))
mgcfp_g = re.findall("G", str(gc_forward_primer))
mgcfp_a = re.findall("A", str(gc_forward_primer))
mgcfp_t = re.findall("T", str(gc_forward_primer))

mrp_c = re.findall("C", str(reverse_primer))
mrp_g = re.findall("G", str(reverse_primer))
mrp_a = re.findall("A", str(reverse_primer))
mrp_t = re.findall("T", str(reverse_primer))

mgcrp_c = re.findall("C", str(gc_reverse_primer))
mgcrp_g = re.findall("G", str(gc_reverse_primer))
mgcrp_a = re.findall("A", str(gc_reverse_primer))
mgcrp_t = re.findall("T", str(gc_reverse_primer))

#Melting point calculations:
#Smaller than 14 nucleotides = Tm= (wA+xT) * 2 + (yG+zC) * 4
#14 nucleotides or larger = Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)

#GC content % calculation for each primer.
try:
    gcc_forward_primer = (len(mfp_g) + len(mfp_c)) / (len(mfp_a) + len(mfp_t) + len(mfp_g) + len(mfp_c)) * 100
except:
    gcc_forward_primer = 0
try:
    gcc_reverse_primer = (len(mrp_g) + len(mrp_c)) / (len(mrp_a) + len(mrp_t) + len(mrp_g) + len(mrp_c)) * 100
except:
    gcc_reverse_primer = 0
if gcf_end > 0:
    gcc_gc_forward_primer = (len(mgcfp_g) + len(mgcfp_c)) / (len(mgcfp_a) + len(mgcfp_t) + len(mgcfp_g) + len(mgcfp_c)) * 100
else:
    gcc_gc_forward_primer = 0

if gcr_end > 0:
    gcc_gc_reverse_primer = (len(mgcrp_g) + len(mgcrp_c)) / (len(mgcrp_a) + len(mgcrp_t) + len(mgcrp_g) + len(mgcrp_c)) * 100
else:
    gcc_gc_reverse_primer = 0

m_gc_primers = 0
if mt.Tm_Wallace(Seq(gc_forward_primer)) or mt.Tm_Wallace(Seq(gc_reverse_primer))  == 0:
    m_gc_primers = 0

#Final output for each primer. GC primer will be blank if no G or C was found in specified bounds.
print(res_sitefw)
print(res_siterv)
gc_style_fw = "green"
gc_style_rv = "green"
gc_seq_fw = ""
gc_seq_rv = ""
gc_fw_five_prime = "5'__"
gc_fw_three_prime = "__3'"
gc_rv_five_prime = "5'__"
gc_rv_three_prime = "__3'"
gc_res_sitefw = ""
gc_res_siterv = ""
if len(gc_forward_primer) < 1:
    gc_style_fw = "red"
    gc_seq_fw = "N/A"
    gc_fw_five_prime = ""
    gc_fw_three_prime = ""
    gc_res_sitefw = "N/A"
else:
    gc_style_fw = gc_forward_primer
    gc_res_sitefw = res_sitefw

if len(gc_reverse_primer) < 1:
    gc_style_rv = "red"
    gc_seq_rv = "N/A"
    gc_rv_five_prime = ""
    gc_rv_three_prime = ""
    gc_res_siterv = "N/A"
else:
    gc_seq_rv = gc_reverse_primer
    gc_res_siterv = res_siterv
    

table = Table(title="[bold]Generated Primers[/bold]", caption=f'[bold]Reference Sequence[/bold]:\n 5\'__{str_seq}__3\'', show_lines=1)
table.add_column("Primer Orientation")
table.add_column("Sequence (5'->3')", overflow="fold")
table.add_column("Restriction Site")
table.add_column("Length (bp)")
table.add_column("Approx. Melting Point (\xb0C)")
table.add_column("GC Content (%)")
table.add_row("[yellow]FORWARD[/yellow]", f"'5__{str(forward_primer)}__3'", res_sitefw, str(len(forward_primer)), str(mt.Tm_Wallace(Seq(forward_primer))), str(round(gcc_forward_primer)), style="green")
table.add_row("[yellow]REVERSE[/yellow]", f"'5__{str(reverse_primer)}__3'", res_siterv, str(len(reverse_primer)), str(mt.Tm_Wallace(Seq(reverse_primer))), str(round(gcc_reverse_primer)), style="green")
table.add_row("[yellow]GC Clamp FORWARD (if applicable)[/yellow]", f"{gc_fw_five_prime}{gc_seq_fw}{gc_fw_three_prime}", gc_res_sitefw, str(len(gc_forward_primer)), str(mt.Tm_Wallace(Seq(gc_forward_primer))), str(round(gcc_gc_forward_primer)), style=gc_style_fw)
table.add_row("[yellow]GC Clamp REVERSE(if applicable)[/yellow]", f"{gc_rv_five_prime}{gc_seq_rv}{gc_rv_three_prime}", gc_res_siterv, str(len(gc_reverse_primer)), str(mt.Tm_Wallace(Seq(gc_reverse_primer))), str(round(gcc_gc_reverse_primer)), style=gc_style_rv)


console = Console()
print('\n\n' + '-' * 120)
console.print(table)
print('\n' + '-' * 120)

quit()