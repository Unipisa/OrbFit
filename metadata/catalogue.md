## DIUNIPI SWH DEPOSITORY/JOURNAL.md Template

Here we propose a template for CATALOGUE<span>.md for Depository of DIUNIPI SWH Acquisition Process.

Each Item in the MATERIAL folder of Depository should have a corresponding record on the CATALOGUE<span>.md with the structure below..

Please note that:
* Name and Surname of actors should be linked to their paragraph in [ACTORS.md](./ACTORS.md) file;
* Items should be linked to the file [inside the repository](./MATERIAL/);
* On the [second part of the Catalogue](./CATALOGUE.md#SW_NAME-Depository-Catalougue-Tree) should be copied the result of the command `tree -a` ;
* Notes are optional;
* *Warehouse:* is optional - should be used only when a physical warehouse is used to store material taken from the *origin*.

Example of Actor link:
~~~
[Name Surname](./ACTORS.md#name-surname)
~~~
Example of Item link:
~~~
[Item Name](./MATERIAL/example_file.zip)
~~~


# SW_NAME Depository Catalougue


* **[Item Name](./MATERIAL/example_file.zip)**
  * *Origin:* 
  * *Warehouse:*
  * *Authors:* [Name Surname](./ACTORS.md#name-surname)
  * *Collectors:* [Name Surname](./ACTORS.md#name-surname)
  * *Description:* 
  * *Notes:*
  
.

# SW_NAME Depository Catalougue Tree


result of `tree -a`  on MATERIAL directory.
