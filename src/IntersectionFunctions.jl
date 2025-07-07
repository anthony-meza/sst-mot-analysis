function return_intersection(x::BoundaryCondition, y::BoundaryCondition)
    x_copy, y_copy = 1. * x, 1. * y
    xnan = isnan.(x_copy.tracer)
    ynan = isnan.(y_copy.tracer)
    x_copy.tracer[ xnan .| ynan] .= NaN
    y_copy.tracer[ xnan .| ynan] .= NaN
    return x, y
end
function plot_map_and_zonal(lon, lat, map_data, zonal_vals;
                        map_norm, cmap, map_title, cb_label,
                        zonal_xticks, zonal_title, savepath, plot_zonal = false)
    fig = figure(figsize=(15,5))
    ax  = fig.add_subplot(131,  projection=ccrs.PlateCarree(central_longitude = 180))

    # — map panel —
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor="#B0B0B0", facecolor = "#B0B0B0")
    ax.coastlines(color="#B0B0B0", alpha=0.3); ax.gridlines(alpha=0.0)

    ax.set_xticks([-180,-90,0,90,180]); ax.set_yticks(-90:15:90)
    mesh = ax.pcolormesh(lon, lat, map_data;
                        cmap=cmap, norm=map_norm,
                        transform=ccrs.PlateCarree())
    ax.set_title(map_title)
    cb = fig.colorbar(mesh, ax=ax;
                      orientation="horizontal",
                      extend="both", fraction=0.03, pad=0.1)
    # cb.set_label(cb_label)
    cb.ax.tick_params(rotation=45)

    bbox = cb.ax.get_position()  # Bounding box in figure coordinates
    x_pos = bbox.x1 + 0.75  # Slightly to the right of the colorbar
    y_pos = (bbox.y0 + bbox.y1) / 2  # Vertically centered
    cb.ax.text(
        x_pos, y_pos, cb_label, 
        transform=cb.ax.transAxes, rotation=0)
    fig.tight_layout()
    
    if plot_zonal
        # # — zonal panel —
        axz = fig.add_subplot(132)
        
        pos1 = ax.get_position()  # returns the bounding box position of the first plot
        height = pos1.height  # Get the height of the first subplot
        
        pos2 = axz.get_position()
        axz.set_position([pos2.x0, pos1.y0, pos2.width /3, height])  # Align the height
    
        
        axz.plot(zonal_vals[:], lat, "o-k"; markersize=5, alpha=0.7)
        axz.set_yticks(-90:15:90); axz.set_ylim(-90, 90)
        axz.tick_params(labelleft=false)    
        axz.set_xticks(zonal_xticks)
        axz.grid(alpha=0.3)
        axz.set_title(zonal_title); axz.set_xlabel(cb_label)
    end
    
    fig.savefig(savepath; dpi=200, bbox_inches="tight")
    if plot_zonal
        return fig, ax, axz
    else 
        return fig, ax
    end
end
# # — Case 1 —


# — Case 2 —

