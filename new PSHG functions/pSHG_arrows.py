import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl


def pSHGArrows(I2, Phi2, allSum, mask, data_path, arrowColourMax):

    # ------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------
    n = 12                  # spacing between arrows (pixels)
    arrowLength = 5       # arrow half-length in pixels
    angularOffset = 0       # degrees

    # ------------------------------------------------------------
    # Colour map for arrow intensity
    # ------------------------------------------------------------
    cmap = plt.cm.jet
    cNorm = colors.Normalize(vmin=0, vmax=arrowColourMax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

    # ------------------------------------------------------------
    # Grayscale image contrast
    # ------------------------------------------------------------
    validImage = allSum[mask != 0]
    
    vmin = np.percentile(validImage, 1)
    vmax = np.percentile(validImage, 99.5)
    
    # Make masked/thresholded regions black
    gray_cmap = plt.cm.gray.copy()
    gray_cmap.set_bad('black')
    
    # ------------------------------------------------------------
    # Set up figure
    # ------------------------------------------------------------
    fig = plt.figure(figsize=(12, 10))
    
    ax = fig.add_axes([0.08, 0.08, 0.75, 0.84])
    axc = fig.add_axes([0.86, 0.08, 0.04, 0.84])
    
    image = np.ma.masked_where(mask == 0, allSum)
    
    ax.imshow(
        image,
        cmap=gray_cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation='nearest'
    )

    # ------------------------------------------------------------
    # Sample arrows
    # ------------------------------------------------------------
    
    
    # Centre of each n x n sampling region
    x = np.arange(n/2 - 0.5, Phi2.shape[1], n)
    y = np.arange(n/2 - 0.5, Phi2.shape[0], n)
    
    for pt0 in x:
        for pt1 in y:
    
            xi = int(round(pt0))
            yi = int(round(pt1))
    
            # Check image boundaries
            if xi >= Phi2.shape[1] or yi >= Phi2.shape[0]:
                continue
    
            # Check mask
            if mask[yi, xi] == 0:
                continue
    
            phi = Phi2[yi, xi]
            intensity = I2[yi, xi]
    
            arrowX = arrowLength * np.cos(
                np.deg2rad(phi + angularOffset)
            )
    
            arrowY = arrowLength * np.sin(
                np.deg2rad(phi + angularOffset)
            )
    
            colorVal = scalarMap.to_rgba(intensity)
    
            ax.arrow(
                pt0, pt1,
                arrowX, arrowY,
                color=colorVal,
                linewidth=2,
                alpha=0.6,
                head_width=1,
                head_length=1,
                length_includes_head=True
            )
    
            ax.arrow(
                pt0, pt1,
                -arrowX, -arrowY,
                color=colorVal,
                linewidth=2,
                alpha=0.6,
                head_width=1,
                head_length=1,
                length_includes_head=True
            )

    # ------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------
    ax.set_xlim(0, Phi2.shape[1])
    ax.set_ylim(Phi2.shape[0], 0)
    ax.axis('off')

    mpl.colorbar.ColorbarBase(
        axc,
        cmap=cmap,
        norm=cNorm,
        orientation='vertical'
    )

    fig.suptitle('PSHG - Arrow colour showing local I2 value')

    plt.savefig(
        data_path + '\\pSHG Arrows.png',
        dpi=400,
        bbox_inches='tight',
        pad_inches=0
    )

    plt.show()