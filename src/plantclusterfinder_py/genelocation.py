import polars as pl
from collections import Counter
from dataclasses import dataclass, field
from plantclusterfinder_py.cluster_datatypes import Conversion

@dataclass
class GeneLoc:
    file_path: str
    data: pl.DataFrame = field(init=False)

    def __post_init__(self):
        """Load data from the specified file path into the Polars DataFrame."""
        self.data = pl.read_csv(self.file_path)  # You can change this to read other formats

    def __len__(self):
        """Return the length of the DataFrame."""
        return len(self.data) # or self.data.height using polars

    def sort_genes(self, chromosome_col: str, start_col: str, end_col: str, strand_col: str):
        """Sort the DataFrame by chromosome name and then conditionally by start or end based on strand."""
        self.data = self.data.sort(
            by=[
                chromosome_col,
                pl.when(pl.col(strand_col) == '+').then(pl.col(start_col)).otherwise(pl.col(end_col))
            ]
        )
    
    def map_to(self, conversion: Conversion, mapping_type: str):
        """Map gene names to conversion.MapT and add a column 'Map_to_T'.
        Args:
            conversion (Conversion): An instance of Conversion containing gene sets (MapG, MapP, MapT).
            mapping_type (str): The type of mapping to use ('G', 'P', or 'T').

        Raises:
            ValueError: If the selected gene set is empty or if mapping_type is invalid.
        """

        if mapping_type == 'G':
            mapping_genes = conversion.MapG
        elif mapping_type == 'P':
            mapping_genes = conversion.MapP
        elif mapping_type == 'T':
            mapping_genes = conversion.MapT
        else:
            raise ValueError("Invalid mapping type! Choose 'G' for Gene, 'P' for Protein, or 'T' for Transcript.")
        
        # Check if the selected gene set is valid
        if not mapping_genes:
            raise ValueError(f"Selected mapping type '{mapping_type}' has no genes to map.")

        # Count occurrences of each gene in the mapping_genes
        gene_count = Counter(mapping_genes)

        # Create a new column with the count of each gene
        column_name = f"Map_to_{mapping_type}"  # Set the column name based on the mapping type
        self.data = self.data.with_columns(
            pl.col("gene name").map_dict(gene_count).fillna(0).cast(pl.Int32).alias(column_name)
        )

    def __repr__(self):
        """Return a string representation of the DataFrame."""
        return self.data.__repr__()



# Example usage
if __name__ == "__main__":
    
    gene_loc = GeneLoc(file_path='file_path')

    # Get the length of the DataFrame
    print(f"Number of entries: {len(gene_loc)}")

    # Sort the genes
    gene_loc.sort_genes(chromosome_col='chromosome', 
                        start_col='gene_start', 
                        end_col='gene_end', 
                        strand_col='strand')
    
    # Display the sorted DataFrame
    print(gene_loc)
