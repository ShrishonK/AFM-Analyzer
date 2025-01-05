#AFM-Analyzer

##Folder Structure
Make sure you create a folder and 3 subfolders with the `.jpk-force-map` files within each file

###Sample Structure:
`main/
├── folder1/           # Root layout with tabs
│   ├── sample1.jpk     # JPK Files
├── folder2/            # Home screen
│   ├── sample2.jpk     # JPK Files
├── folder3/     # Theme context
│   ├── sample3.jpk     # JPK Files`

##Steps to Run the Python Program
1.Use the CD function to enter the directory of the outer folder.

2.To install the required packages, run "pip install -r requirements.txt"

3.To run export the CSV, enter the full file path and then comment out the line of code that shows the plot.

##Files Returned:
- .csv file with the Filename, Contact Point, Turnaround Point, and Area
- Graphs for each curve
- Average Plot Graph for each folder 
