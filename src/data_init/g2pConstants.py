"CONSTANTS"

import os

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
data_dir = os.path.join(src, os.pardir, "data")
archive_dir = os.path.join(data_dir, "archive")
src, data_dir = map(os.path.abspath, (src, data_dir))
res_dict = {
	"R": "RESISTANT", "R*": "LIKELY RESISTANT", "R?": "POSSIBLY RESISTANT",
	"S": "SENSITIVE", "S*": "LIKELY SENSITIVE", "?": "AMBIGUOUS",
	"H": "HYPER-SENSITIVE", "H*": "LIKELY HYPER-SENSITIVE"
}
