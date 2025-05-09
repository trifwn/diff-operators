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
    plot_results("test_case1.dat")