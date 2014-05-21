Perception of Statistical Graphics
========================================================
author: Susan VanderPlas
date: May 27, 2014
css: Prelim.css





Outline
========================================================
type:sub-section
- Importance of understanding perception
- Sine Illusion
- Visual Aptitude and Graphical Inference
- Salience of Graph Features

Signs of the Sine Illusion
========================================================
type:sub-section
## Why we need to care
<img src="Prelim-figure/original.png" title="plot of chunk original" alt="plot of chunk original" style="display: block; margin: auto;" />

Example of the Sine Illusion
========================================================
8-hour Average Ozone Levels in Houston, TX by temperature at Hobby Airport

<img src="Prelim-figure/example-sineillusion.png" title="plot of chunk example-sineillusion" alt="plot of chunk example-sineillusion" style="display: block; margin: auto;" />

Explaining the Illusion
========================================================
<img src="Prelim-figure/ribbon-illusion.png" title="Perspective plot similar to the sine illusion" alt="Perspective plot similar to the sine illusion" style="display: block; margin: auto;" />
Perspective plot of a three-dimensional image similar to the sine illusion


Explaining the Illusion
========================================================
<img src="Prelim-figure/ribbon-illusion2.png" title="plot of chunk ribbon-illusion2" alt="plot of chunk ribbon-illusion2" style="display: block; margin: auto;" />
Perspective plot of the same data, with a vanishing point closer to infinity

Explaining the Illusion
========================================================
<img src="Prelim-figure/originalgrid.png" title="plot of chunk originalgrid" alt="plot of chunk originalgrid" style="display: block; margin: auto;" />
The illusion doesn't disappear with grid lines, but does disappear when the context is removed

Geometry of the Sine Illusion
========================================================
<img src="Prelim-figure/illusion-geometry.png" title="plot of chunk illusion-geometry" alt="plot of chunk illusion-geometry" width="90%" height="auto" style="display: block; margin: auto;" />
 
Geometry of the Sine Illusion
========================================================
- We perceive the orthogonal width of the implied surface
- The orthogonal width is a function of the x and y range as well as the aspect ratio of the plot. 
- The perceived orthogonal width is also a function of the slope of the line tangent to the underlying function curve. 

Correcting the Illusion
=======================================================
- Remove the underlying function (plot the curve and the residuals separately)
- Change the plotted line length (or spread) so that the **perceived** orthogonal width corresponds to the **original** (data) line length
- Reparameterize the x-axis in terms of the slope, so that the absolute slope doesn't change

Correcting the Illusion 
=======================================================
## Trend Removal
<img src="Prelim-figure/cleveland.png" title="plot of chunk cleveland" alt="plot of chunk cleveland" style="display: block; margin: auto;" />

Correcting the Illusion
=======================================================
## Trend Removal
<img src="Prelim-figure/cleveland2.png" title="plot of chunk cleveland2" alt="plot of chunk cleveland2" style="display: block; margin: auto;" />

Correcting the Illusion
=======================================================
## X-axis Transformation
<img src="Prelim-figure/xaxisdemo.png" title="plot of chunk xaxisdemo" alt="plot of chunk xaxisdemo" style="display: block; margin: auto;" />

Correcting the Illusion
=======================================================
## X-axis Transformation
<img src="Prelim-figure/xaxisdemoweight.png" title="plot of chunk xaxisdemoweight" alt="plot of chunk xaxisdemoweight" style="display: block; margin: auto;" />

Correcting the Illusion
=======================================================
## Y-axis Transformation
<img src="Prelim-figure/y-generalcorrectioncartoon.png" title="plot of chunk y-generalcorrectioncartoon" alt="plot of chunk y-generalcorrectioncartoon" style="display: block; margin: auto;" />
Correcting the Illusion
=======================================================
## Y-axis Transformation
<img src="Prelim-figure/y-linearcorrectioncartoon.png" title="plot of chunk y-linearcorrectioncartoon" alt="plot of chunk y-linearcorrectioncartoon" style="display: block; margin: auto;" />

Correcting the Illusion
=======================================================
## Y-axis Transformation
<img src="Prelim-figure/ycorrection.png" title="plot of chunk ycorrection" alt="plot of chunk ycorrection" width="100%" style="display: block; margin: auto;" />
