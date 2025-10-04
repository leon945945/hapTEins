**software**
- pipeline-detectTE.sh:

	line 10, 13, 14, 17

- getPos.py:

	line 44, 65

- fltDiscordant.py:

	line 55, 90, 105, 106, 112, 113, 114

**paramater**
- getPos.py:
	line 10: TE names appended in genome file

	line 12: TE length get from TE name

	line 56: softclip distance

	line 67: read depth among insertion site, 10-100

- fltDiscordant.py:
	line 26: keep sites with depth 10-50
  
	line 68: softclip distance < 20 && reads number > 10

	line 78: depth of candidate insertion site ranging from 10-50
