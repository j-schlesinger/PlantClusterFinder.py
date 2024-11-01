import polars as pl
from dataclasses import dataclass, field

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
