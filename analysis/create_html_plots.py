#!/usr/bin/env python3

import json
import html

# Read the analysis data
with open('peak_analysis_data.json', 'r') as f:
    data = json.load(f)

# Create HTML visualization
html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GSE130638 ATAC-seq Analysis Results</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        h1, h2 { color: #2c3e50; text-align: center; }
        .plot-container { margin: 30px 0; padding: 20px; border: 1px solid #ddd; border-radius: 5px; background-color: #fafafa; }
        .stats-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 12px; text-align: center; }
        .stats-table th { background-color: #3498db; color: white; }
        .stats-table tr:nth-child(even) { background-color: #f2f2f2; }
        .summary { background-color: #e8f6f3; padding: 15px; border-radius: 5px; margin: 20px 0; }
        .highlight { color: #e74c3c; font-weight: bold; }
    </style>
</head>
<body>
    <div class="container">
        <h1>GSE130638 ATAC-seq Analysis Results</h1>
        <h2>Comparing G82 (Parental) vs G82R (Resistant) Brain Tumor Cells</h2>
        
        <div class="summary">
            <h3>Key Findings Summary</h3>
            <ul>
                <li>G82R shows <span class="highlight">10.9% more peaks</span> than G82 parental cells</li>
                <li>Mean enrichment score increased by <span class="highlight">14.9%</span> in resistant cells</li>
                <li>Mean signal intensity increased by <span class="highlight">18.9%</span> in resistant cells</li>
                <li>Peak widths are slightly larger in G82R (4.3% increase)</li>
                <li>Chromosome 3, 6, 9, and 17 show the largest increases in accessible regions</li>
            </ul>
        </div>

        <!-- Peak Count Comparison -->
        <div class="plot-container">
            <h3>1. Peak Count Comparison</h3>
            <div id="peakCountPlot"></div>
        </div>

        <!-- Peak Width Statistics -->
        <div class="plot-container">
            <h3>2. Peak Width Statistics</h3>
            <table class="stats-table">
                <tr>
                    <th>Statistic</th>
                    <th>G82 (Parental)</th>
                    <th>G82R (Resistant)</th>
                    <th>Change</th>
                </tr>"""

# Add width statistics table
width_stats = data['width_stats']
for stat in ['mean', 'median', 'q1', 'q3']:
    g82_val = width_stats['G82'][stat]
    g82r_val = width_stats['G82R'][stat]
    change = ((g82r_val - g82_val) / g82_val) * 100
    html_content += f"""
                <tr>
                    <td>{stat.upper()}</td>
                    <td>{g82_val:.1f} bp</td>
                    <td>{g82r_val:.1f} bp</td>
                    <td>{change:+.1f}%</td>
                </tr>"""

html_content += """
            </table>
        </div>

        <!-- Chromosome Distribution -->
        <div class="plot-container">
            <h3>3. Chromosome Distribution Changes</h3>
            <div id="chrDistPlot"></div>
        </div>

        <!-- Enrichment Score Comparison -->
        <div class="plot-container">
            <h3>4. Enrichment Score Statistics</h3>
            <table class="stats-table">
                <tr>
                    <th>Statistic</th>
                    <th>G82 (Parental)</th>
                    <th>G82R (Resistant)</th>
                    <th>Change</th>
                </tr>"""

# Add enrichment statistics table
enrich_stats = data['enrichment_stats']
for stat in ['mean', 'median', 'q1', 'q3']:
    g82_val = enrich_stats['G82'][stat]
    g82r_val = enrich_stats['G82R'][stat]
    change = ((g82r_val - g82_val) / g82_val) * 100
    html_content += f"""
                <tr>
                    <td>{stat.upper()}</td>
                    <td>{g82_val:.2f}</td>
                    <td>{g82r_val:.2f}</td>
                    <td>{change:+.1f}%</td>
                </tr>"""

html_content += """
            </table>
        </div>

        <!-- Signal Intensity Comparison -->
        <div class="plot-container">
            <h3>5. Signal Intensity Statistics</h3>
            <table class="stats-table">
                <tr>
                    <th>Statistic</th>
                    <th>G82 (Parental)</th>
                    <th>G82R (Resistant)</th>
                    <th>Change</th>
                </tr>"""

# Add signal statistics table
signal_stats = data['signal_stats']
for stat in ['mean', 'median', 'q1', 'q3']:
    g82_val = signal_stats['G82'][stat]
    g82r_val = signal_stats['G82R'][stat]
    change = ((g82r_val - g82_val) / g82_val) * 100
    html_content += f"""
                <tr>
                    <td>{stat.upper()}</td>
                    <td>{g82_val:.1f}</td>
                    <td>{g82r_val:.1f}</td>
                    <td>{change:+.1f}%</td>
                </tr>"""

html_content += """
            </table>
        </div>

        <!-- Peak Overlap Analysis -->
        <div class="plot-container">
            <h3>6. Peak Overlap Analysis</h3>
            <div id="overlapPlot"></div>
            <div class="summary">
                <p><strong>Note:</strong> The very low overlap (0.1%) suggests this analysis used exact coordinate matching. 
                In practice, ATAC-seq peaks should be analyzed with proper genomic overlap tools that account for 
                nearby or partially overlapping regions.</p>
            </div>
        </div>

        <div class="summary">
            <h3>Biological Interpretation</h3>
            <p>The G82R resistant cells show a pattern of <strong>increased chromatin accessibility</strong> compared to parental G82 cells. 
            This suggests that resistance to CAY10566 may involve epigenetic reprogramming that opens new regulatory regions, 
            potentially enabling alternative gene expression programs. The increased signal intensities and peak counts indicate 
            a more "open" chromatin state that could facilitate transcriptional plasticity during resistance development.</p>
        </div>
    </div>

    <script>
        // Peak Count Comparison
        var peakCountData = [{
            x: ['G82 (Parental)', 'G82R (Resistant)'],
            y: [""" + str(data['peak_counts']['G82']) + """, """ + str(data['peak_counts']['G82R']) + """],
            type: 'bar',
            marker: { color: ['#3498db', '#e74c3c'] },
            text: [""" + str(data['peak_counts']['G82']) + """, """ + str(data['peak_counts']['G82R']) + """],
            textposition: 'outside'
        }];
        
        var peakCountLayout = {
            title: 'Total ATAC-seq Peaks',
            yaxis: { title: 'Number of Peaks' },
            showlegend: false
        };
        
        Plotly.newPlot('peakCountPlot', peakCountData, peakCountLayout);

        // Chromosome Distribution
        var chrNames = [""" + ', '.join([f"'{item['chr']}'" for item in data['chr_distribution']]) + """];
        var pctChanges = [""" + ', '.join([str(item['pct_change']) for item in data['chr_distribution']]) + """];
        var colors = pctChanges.map(x => x < 0 ? '#e74c3c' : '#27ae60');
        
        var chrData = [{
            x: pctChanges,
            y: chrNames,
            type: 'bar',
            orientation: 'h',
            marker: { color: colors },
            text: pctChanges.map(x => x.toFixed(1) + '%'),
            textposition: 'outside'
        }];
        
        var chrLayout = {
            title: 'Percentage Change in Peak Counts by Chromosome',
            xaxis: { title: 'Percentage Change (%)' },
            yaxis: { title: 'Chromosome' },
            showlegend: false
        };
        
        Plotly.newPlot('chrDistPlot', chrData, chrLayout);

        // Peak Overlap Visualization
        var overlapData = [{
            labels: ['G82 Specific', 'G82R Specific', 'Overlapping'],
            values: [""" + str(data['overlap']['g82_unique']) + """, """ + str(data['overlap']['g82r_unique']) + """, """ + str(data['overlap']['total']) + """],
            type: 'pie',
            marker: { colors: ['#3498db', '#e74c3c', '#27ae60'] }
        }];
        
        var overlapLayout = {
            title: 'Peak Overlap Distribution'
        };
        
        Plotly.newPlot('overlapPlot', overlapData, overlapLayout);
    </script>
</body>
</html>
"""

# Write HTML file
with open('ATAC_seq_analysis_report.html', 'w') as f:
    f.write(html_content)

print("Interactive HTML report created: ATAC_seq_analysis_report.html")
print("Open this file in a web browser to view the interactive plots and analysis results.")