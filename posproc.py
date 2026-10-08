import numpy as np
import matplotlib.pyplot as plt

def main():
    # Load data from the tracking file
    # Columns: TIME_VAL, disp_x, disp_y, trac_x, trac_y
    data = np.loadtxt('disp_trac.dat')
    
    if data.size == 0:
        print("Error: disp_trac.dat is empty.")
        return
        
    # Depending on whether there is one row or multiple, data might be 1D or 2D
    if data.ndim == 1:
        data = data.reshape(1, -1)
        
    time_vals = data[:, 0]
    disp_y = data[:, 2]
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_vals, disp_y, marker='o', linestyle='-', color='b', label='Deslocamento em Y')
    
    # Add labels and title
    plt.xlabel('Tempo (s)', fontsize=12)
    plt.ylabel('Deslocamento Y', fontsize=12)
    plt.title('Deslocamento na direção Y ao longo do tempo', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Save the figure to a PNG file
    output_filename = 'deslocamento_y_tempo.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Figura salva com sucesso em: {output_filename}")

if __name__ == '__main__':
    main()
