// Copied from: https://git.sr.ht/~freestatelabs/csv/tree/main/item/csv.h

#ifndef CSV_H 
#define CSV_H 

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define MAX_LINE_LENGTH 1000


// Data is stored in row-major format
struct CsvData {
    int nrows, ncols; 
    double *values;
};


// Convenience function for freeing the data structure
// (not needed)
static void 
free_csv_data(struct CsvData data) {
    free(data.values);
}


// Prints the data to the terminal
static void 
print_csv_data(struct CsvData data) {
    for (int i=0; i<data.nrows; i++) {
        for (int j=0; j<data.ncols; j++) {
            printf("%.6f ", data.values[data.ncols*i + j]);
        }
        printf("\n");
    }
}

// Maintains list of valid characters for data
static bool 
valid_char(char c) {
    char valid[15] = "0123456789eE.-";
    for (int i=0; i<15; i++) {
        if (c==valid[i]) return true; 
    }
    return false;
}

// Parse a line from the file and store in line_data
// Manually checks every character (slower)
static void 
parse_line_manual(const char s[], double *line_data, int ncols, char delimiter) { 

    int col = 0; 
    char word[MAX_LINE_LENGTH]; 
    int word_index = 0; 
    int line_index = 0; 

    for (int i=0; s[i] != '\0'; i++) {
        char c = s[i];
        if (valid_char(c)) {
            word[word_index] = c; 
            word_index++;
        }
        if (c == delimiter) {
            word[word_index] = '\0';
            line_data[line_index] = atof(word); 
            line_index++;
            word_index = 0; 
            col++; 
        }
        if (ncols == col) break;
    }
    // Last value might not have a delimiter after it
    word[word_index] = '\0';
    line_data[line_index] = atof(word);
}


// Parse a single line from the file and store in `line_data`
inline static void 
parse_line(char s[], double *line_data, int ncols, char delimiter) {
    int col = 0; 
    int i = 0; 
    char *start_ptr = s;
    char *end_ptr;
    while (col < ncols) {
        if (s[i] == delimiter || s[i] == '\0') {
            end_ptr = s + i;
            line_data[col] = strtod(start_ptr, &end_ptr);
            start_ptr = end_ptr+1;
            col++;
        }
        if (s[i] == '\0') break;
        else i++;
    }
}


// Read a csv file, allocate memory, and parse the data from the file
static struct CsvData 
read_csv(
    const char fn[], int ncols, int nrows, int skiprows, char delimiter
) {

    struct CsvData data = {
        .nrows=nrows, 
        .ncols=ncols, 
        .values=(double*)calloc(nrows*ncols, sizeof(double))
    }; 

    char s[MAX_LINE_LENGTH];
    double *line_data = (double*)calloc(ncols, sizeof(double));

    int row = 0;            // of data 
    int line = 0;           // of the file; different depending on how many rows are skipped

    FILE *file = fopen(fn, "r");
    if (file == NULL) return (struct CsvData){.nrows=0, .ncols=0, .values=NULL};
    while (fgets(s, sizeof(s), file)) {
        if (line >= skiprows) {
            parse_line(s, line_data, ncols, delimiter);
            for (int i=0; i<ncols; i++) {
                data.values[ncols*row+i] = line_data[i];
            }
            row++;   
        }
        line++;

        // Reallocate if the file is larger than expected
        if (row >= nrows) {
            nrows *= 2; 
            double *new_ptr = (double*)realloc(data.values, nrows*ncols*sizeof(double));
            if (new_ptr != NULL) {
                data.values = new_ptr;
                data.nrows = nrows;
            }
            else {
                printf("ERROR! Could not realloc() memory in read_csv()!\n");
            }
        }

    }

    free(line_data);
    fclose(file); 
    nrows = row;
    double *new_ptr = (double*)realloc(data.values, nrows*ncols*sizeof(double));
    if (new_ptr != NULL) {
        data.values = new_ptr;
        data.nrows = nrows;
    }

    return data;
}


// Write a matrix out to a csv file
static void 
write_csv(
    char *fn, struct CsvData data, char delimiter, char header[]
) {

    FILE *file = fopen(fn, "w"); 
    if (file == NULL) return;
    fprintf(file, header);
    fprintf(file, "\n");
    for (int i=0; i<data.nrows; i++) {
        for (int j=0; j<data.ncols; j++) {
            fprintf(file, "%.17g", data.values[i*data.ncols+j]);
            fprintf(file,",");
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

#endif

