#!/usr/bin/env python

import pandas as pd
from datetime import datetime

daily_filename = 'runoff._.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.csv'
daily_filename = 'smb_rec._.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.csv'
monthly_filename = daily_filename.replace('.DD.','.MM.')

df = pd.read_csv(daily_filename, header=0, parse_dates=['#date'], index_col='#date')
monthly_sum = df.resample('M').sum()
monthly_sum.to_csv(monthly_filename)

#import csv
#dts = list()
#values = list()
#with open('runoff._.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.csv') as csv_file:
#   csv_reader = csv.reader(csv_file, delimiter=',')
#   next(csv_reader, None)
#   for irow, row in enumerate(csv_reader):
#      dts.append(datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S'))
#      values.append(float(row[1]))

