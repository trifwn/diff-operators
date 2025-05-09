import numpy as np
import matplotlib.pyplot as plt
import os

def plot_results(filename):
    if not os.path.exists(filename):
        print(f"File {filename} not found.")
        return
    
    # Load data
    data = np.loadtxt(filename, skiprows=1)
    x = data[:, 0]
    f = data[:, 1]
    df_exact = data[:, 2]
    df_compact = data[:, 3]
    df_fd = data[:, 4]
    error_compact = data[:, 5]
    error_fd = data[:, 6]
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot derivatives
    ax1.plot(x, df_exact, 'k-', linewidth=2, label='Exact')
    ax1.plot(x, df_compact, 'r--', linewidth=1.5, label='Compact')
    ax1.plot(x, df_fd, 'b:', linewidth=1.5, label='Finite Difference')
    ax1.set_title('Derivative Comparison')
    ax1.set_xlabel('x')
    ax1.set_ylabel('df/dx')
    ax1.legend()
    ax1.grid(True)
    
    # Plot errors on logarithmic scale
    ax2.semilogy(x, error_compact, 'r-', linewidth=1.5, label='Compact Error')
    ax2.semilogy(x, error_fd, 'b-', linewidth=1.5, label='FD Error')
    ax2.set_title('Error Comparison (Log Scale)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('Absolute Error')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig(f"{os.path.splitext(filename)[0]}_plot.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    # Test case
    import os
    test_case_dir = os.path.dirname(__file__)
    # Find all test_case.*.dat files in the directory
    test_case_files = [f for f in os.listdir(test_case_dir) if f.startswith("test_case") and f.endswith(".dat")]
    if not test_case_files:
        print("No test_case.*.dat files found.")
    else:
        for test_case_file in test_case_files:
            print(f"Plotting results for {test_case_file}...")
            plot_results(test_case_file)