## Usage 

```
Usage: words.py [OPTIONS] QUERY_SEQUENCE DATA_SEQUENCE

Options:
  -t, --threshold INTEGER  Threshold score for k-mer pair match.
  -k INTEGER               k-mer length
  --help                   Show this message and exit.
```

```
$ python words.py ACEDECADE REDACEACECEDKL -k 3 -t 10 | sort | uniq -c
      2 ACE ACE 18
      1 ACE ECE 13
      1 CAD CED 14
      1 CED CEA 12
      1 CED CEC 11
      1 CED CED 20
      1 DEC CEC 11
      1 DEC DAC 14
      1 ECA ECE 13
      1 EDE EDK 12
```
