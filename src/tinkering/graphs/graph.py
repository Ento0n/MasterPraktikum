import plotly.graph_objects as go

# Define the x values and original y values
x_values = ["UniProt", "GenBank accession", "GenBank CDSs", "Correct reconstructed"]
# original_y_values = [15224, 15201, 15155, 12321]
original_y_values = [15224, 15201, 15155, 14469]

# Convert the y values to percentages
y_values = [value / original_y_values[0] * 100 for value in original_y_values]

# Create a line graph
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=x_values,
    y=y_values,
    mode='lines+markers+text',
    text=[f"{value:.2f}%" for value in y_values],
    textposition='top center',
    name='Values',
    line=dict(color='blue'),
    marker=dict(size=8)
))

fig.update_yaxes(rangemode="tozero")

# Update the layout to remove y-axis values, title, and x-axis description
fig.update_layout(
    title='',
    xaxis_title='',
    # yaxis=dict(showticklabels=False),
    xaxis=dict(tickmode='array', tickvals=list(range(len(x_values))), ticktext=x_values),
    template='plotly_white'
)

# Save the figure as an HTML file
fig.write_image("lost_entries_human.jpg")

# Show the figure
fig.show()
