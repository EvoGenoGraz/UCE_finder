
__Testing the script__

```bash
bash$ python ../devel/find_UCEs.py
```

_OUTPUT_
```bash

Find UCE candidates in Mummer output
version: v.1.0


###############
# Disclaimer: #

Tested in Python 2.7.12
Requires Biopython

The script is currently in the beta phase. Any feedback is much appreciated.

##############



usage: find_UCEs.py [-h] [--mummer_raw <FILE>] [--ref_fasta <FILE>]
                    [--que_fasta <FILE>] [--min_length <INT>]
                    [--max_merge_distance <INT>] [--unique_distance <INT>]
                    [--dist_from_end <INT>] [--max_merge_distance_diff <INT>]
                    [--extend <INT>] [--output_length_dists]


```


__Test the usage__

```bash
bash$ python ../devel/find_UCEs.py -h
```

_OUTPUT_
```bash
usage: find_UCEs.py [-h] [--mummer_raw <FILE>] [--ref_fasta <FILE>]
                    [--que_fasta <FILE>] [--min_length <INT>]
                    [--max_merge_distance <INT>] [--unique_distance <INT>]
                    [--dist_from_end <INT>] [--max_merge_distance_diff <INT>]
                    [--extend <INT>] [--output_length_dists]

Find UCE candidates in Mummer output
version: v.1.0

###############
# Disclaimer: #

Tested in Python 2.7.12
Requires Biopython

The script is currently in the beta phase. Any feedback is much appreciated.

##############

optional arguments:
  -h, --help            show this help message and exit
  --mummer_raw <FILE>   Mummer output file.
  --ref_fasta <FILE>    Reference assembly in fasta format.
  --que_fasta <FILE>    Query assembly in fasta format.
  --min_length <INT>    Minimum length of candidates to be retained.
                        Default=0, i.e. no filter
  --max_merge_distance <INT>
                        Maximum allowed distance between Mummer matches on
                        reference sequence to be merged into a single
                        candidate. default=0, i.e. don't merge.
  --unique_distance <INT>
                        Minimum distance between two adjacent Mummer matches
                        on reference sequence to be consider the matches as
                        separate candidates. default=1000
  --dist_from_end <INT>
                        Minimum distance of Mummer match from end of Reference
                        contig/scaffold to be considered as candidate.
                        default=10
  --max_merge_distance_diff <INT>
                        Maximum difference of distances between adjacent
                        matches to be merged on the reference and query.
                        default=0.
  --extend <INT>        Extract candidates with additional up/downstream
                        sequence of this length into separate file. Default=0,
                        i.e. don't extract extended sequences. Default=0
  --output_length_dists
                        Output to screen summaries of the candidate length
                        distributions.

   examples: 
	# Minimum command
	./find_UCEs.py --mummer_raw mummer.out --ref_fasta ref.fa --que_fasta que.fa 
	
	# apply some reasonable filters
	./find_UCEs.py --mummer_raw mummer.out --ref_fasta ref.fa --que_fasta que.fa --min_length 100 --max_merge_distance 40 --unique_distance 1000 --dist_from_end 10 --max_merge_distance_diff 5 --extend 500 --output_length_dists
	

```

__Test run__

```bash
bash$ python ../devel/find_UCEs.py --mummer test_data/LvSf80.out --ref_fasta test_data/Lvar_reduced.fasta --que_fasta test_data/Sfran_reduced.fasta --min_length 100 --max_merge_distance 40 --unique_distance 1000 --dist_from_end 10 --max_merge_distance_diff 5 --extend 500 --output_length_dists
```

_OUTPUT_
```bash
LvSf80.out.tsv
Total number of MUMMER matches: 1482
Filter 'multiple query hits with identical start position on reference' removed 26 Mummer hits
Leaves us with 1456 hits to process

MERGING


Length distributions for candidate groups:

unique	153
80 5
81 6
82 5
83 4
84 6
85 5
86 2
87 3
88 2
89 5
90 4
91 2
92 2
93 5
94 4
95 2
97 4
98 4
99 2
100 4
101 1
102 4
103 2
104 2
105 2
106 3
107 1
108 3
109 2
111 1
112 4
113 2
114 3
115 3
116 1
120 2
121 1
124 1
125 1
127 3
128 2
129 1
131 2
132 6
133 2
134 1
136 1
137 2
140 1
141 1
152 1
154 1
157 1
166 1
170 1
175 1
176 1
180 1
195 1
197 1
209 1
221 1
247 2
309 1
392 1

singlematch	1032
80 73
81 54
82 44
83 38
84 32
85 46
86 43
87 37
88 31
89 35
90 24
91 21
92 25
93 21
94 18
95 37
96 25
97 19
98 20
99 12
100 17
101 24
102 15
103 16
104 12
105 9
106 14
107 16
108 13
109 13
110 9
111 8
112 9
113 12
114 6
115 13
116 7
117 7
118 7
119 10
120 4
121 7
122 5
123 9
124 5
125 8
126 5
127 3
128 8
129 6
130 3
131 1
134 4
135 1
136 3
138 4
139 2
140 1
141 2
142 1
143 1
144 3
145 1
146 6
147 3
148 3
149 2
150 4
151 3
153 2
154 2
155 1
157 1
160 1
161 2
162 3
163 1
164 1
165 1
168 1
169 1
173 1
176 2
177 1
179 1
185 2
190 1
191 1
194 1
195 1
201 1
243 1
274 1
438 1

merged	29
171 1
172 1
179 1
180 1
189 1
190 1
191 1
192 1
194 1
196 2
197 1
207 1
210 1
211 2
212 1
216 1
219 1
220 1
223 1
232 1
233 1
252 1
282 1
326 1
328 1
339 1
535 1

### SUMMARY ###

Candidates dropped:
too-close	7
skip	180
overlap-filter	6
merge-dist-diff	1
position	11
irregular-merge	1

Candidates retained:
unique	153
singlematch	1032
merged	29
_________________
Grand TOTAL: 1214
=================

#### PRODUCING OUTPUT FILES ####

WRITING GLOBAL OUTPUTS with prefix: 'LvSf80.out'
WRITING: Lvar_reduced.candidate_UCEs-minlength100.fasta
WRITING: Sfran_reduced.candidate_UCEs-minlength100.fasta
WRITING: Lvar_reduced.candidate_UCEs-minlength100-extended500.fasta
WRITING: Sfran_reduced.candidate_UCEs-minlength100-extended500.fasta
```
