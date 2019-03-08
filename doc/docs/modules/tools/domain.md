### Domain package

This package allow you to manipulate a domain table

#### Quick Start

This package accepts a domain table, or more broadly a pandas dataframe.

A simple input table consist of:

| ID           |domain  |start   |end|evalue|
|--------------|--------|------- |---|------|
|WP_011523244.1|Reg_prop|270|290|1.00e-02|
|WP_011523244.1|Reg_prop|392|412|5.20e-03|
|WP_011523244.1|Y_Y_Y   |690|752|2.10e-15|
|WP_011523244.1|HisKA_3 |799|864|5.70e-13|
|WP_011523244.1|HATPase_c|915|1004|1.10e-07|

This is a output of rotifer.tools.search().parser().hmmscan_parsed or using architecture2table.

Loading the data:
```python3
from rotifer import domain

# The df is a pandas dataframe containing the domain table
dom = domain(df)
```

After loading the data, the corresponding domain dataframe can be accessed using the domain attribute or using the object name. The main difference is that the domain attribute is a pandas dataframe while the dom (in the example) is a object belonging to the class domain.

```python3
dom = domain(df)
# Similar output

dom.domain
>>>
                  ID                    domain      start  end   evalue
               0  WP_011523244.1         Reg_prop    270   290   1.000000e-02
               1  WP_011523244.1         Reg_prop    392   412   5.200000e-03
               2  WP_011523244.1            Y_Y_Y    690   752   2.100000e-15
               3  WP_011523244.1          HisKA_3    799   864   5.700000e-13
               4  WP_011523244.1        HATPase_c    915  1004   1.100000e-07

dom
>>>
                  ID                    domain      start  end   evalue
               0  WP_011523244.1         Reg_prop    270   290   1.000000e-02
               1  WP_011523244.1         Reg_prop    392   412   5.200000e-03
               2  WP_011523244.1            Y_Y_Y    690   752   2.100000e-15
               3  WP_011523244.1          HisKA_3    799   864   5.700000e-13
               4  WP_011523244.1        HATPase_c    915  1004   1.100000e-07
```

Therefore all pandas dataframe operations are possible when we access the domain attribute.

#### Domain length

To measure the domain length is simple using the domain_len method.

Using the above example we have:
```python3
dom.domain_len()
                  ID                    domain      start  end   evalue     5
               0  WP_011523244.1         Reg_prop    270   290   1.00e-02   21
               1  WP_011523244.1         Reg_prop    392   412   5.20e-03   21
               2  WP_011523244.1            Y_Y_Y    690   752   2.10e-15   63
               3  WP_011523244.1          HisKA_3    799   864   5.70e-13   66
               4  WP_011523244.1        HATPase_c    915  1004   1.10e-07   90
```

Some parameters includes, the start_col and end_col, indicating the domain "start" and "end". The col_name is the ouput column. The inplace parameter if True will modify the domain dataframe.

To change to column name with the calculated size use col_name.
```python3
dom.domain_len(col_name = 'Size')
                  ID                    domain      start  end   evalue     Size
               0  WP_011523244.1         Reg_prop    270   290   1.00e-02   21
               1  WP_011523244.1         Reg_prop    392   412   5.20e-03   21
               2  WP_011523244.1            Y_Y_Y    690   752   2.10e-15   63
               3  WP_011523244.1          HisKA_3    799   864   5.70e-13   66
               4  WP_011523244.1        HATPase_c    915  1004   1.10e-07   90
```

#### Filter by size

col => Size col
max_length = 
min_length
inplace

#### Distribution
Calculate domain distribution
col domain column
col_name 
merge
sort,
ascending
frequency
frequency_col_name

#### plot size

#### percentiles

#### Add sequence

#### slice sequence

seq len


