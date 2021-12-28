
//Time frame                            
//Pre-fire
var prefire_start = '2020-12-20';   
var prefire_end = '2021-07-27';

//Post-fire 
var postfire_start = '2021-08-10'; 
var postfire_end = '2021-12-09';


//Satellite platform (S2/L8 !L8 refresh rate 16 days!)
var platform = 'S2';               


// Satellite platform and dates
if (platform == 'S2' | platform == 's2') {
  var ImCol = 'COPERNICUS/S2';
  var pl = 'Sentinel-2';
} else {
  var ImCol = 'LANDSAT/LC08/C01/T1_SR';
  var pl = 'Landsat 8';
}
print(ee.String('Data selected for analysis: ').cat(pl));
print(ee.String('Fire incident occurred between ').cat(prefire_end).cat(' and ').cat(postfire_start));

// Location (don't forget to import geometry asset/define geometry)
var area = ee.FeatureCollection(geometry);

Map.centerObject(area);


//--
//Sat imagery + collection

var imagery = ee.ImageCollection(ImCol);
var prefireImCol = ee.ImageCollection(imagery
    .filterDate(prefire_start, prefire_end)
    .filterBounds(area));
    
var postfireImCol = ee.ImageCollection(imagery
    .filterDate(postfire_start, postfire_end)
    .filterBounds(area));

print("Pre-fire Image Collection: ", prefireImCol); 
print("Post-fire Image Collection: ", postfireImCol);

//--
//Masks
//Cloud (sentinel + try FAI mask code)
function maskS2sr(image) {
  //10 = clouds, 11 = cirrus 
  var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = ee.Number(2).pow(11).int();
  var qa = image.select('QA60');
  //0 = clear weather
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

//L8 mask cloud from pix quality (for shadows etc)
function maskL8sr(image) {
  //3 = cloud shadows, 5 = clouds
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var snowBitMask = 1 << 4; //not necessary for polbnda_adm_west
 
  var qa = image.select('pixel_qa');
  //0 = clear weather
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
      .and(qa.bitwiseAnd(snowBitMask).eq(0));

  return image.updateMask(mask)
      .select("B[0-9]*")
      .copyProperties(image, ["system:time_start"]);
}

//Platform specific cloud mask
if (platform == 'S2' | platform == 's2') {
  var prefire_CM_ImCol = prefireImCol.map(maskS2sr);
  var postfire_CM_ImCol = postfireImCol.map(maskS2sr);
} else {
  var prefire_CM_ImCol = prefireImCol.map(maskL8sr);
  var postfire_CM_ImCol = postfireImCol.map(maskL8sr);
}
var geometry = geometry;    //modify geometry here

//--
//Mosaic data pre- post-fire
var pre_mos = prefireImCol.mosaic().clip(area);
var post_mos = postfireImCol.mosaic().clip(area);

var pre_cm_mos = prefire_CM_ImCol.mosaic().clip(area);
var post_cm_mos = postfire_CM_ImCol.mosaic().clip(area);

print("Pre-fire True Color Image: ", pre_mos); 
print("Post-fire True Color Image: ", post_mos);

//--
// NBR
// Platform specific NBR = (NIR-SWIR2) / (NIR+SWIR2)
if (platform == 'S2' | platform == 's2') {
  var preNBR = pre_cm_mos.normalizedDifference(['B8', 'B12']);
  var postNBR = post_cm_mos.normalizedDifference(['B8', 'B12']);
} else {
  var preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7']);
  var postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
}

//--
//Difference between pre- post- images
// dNBR
var dNBR_unscaled = preNBR.subtract(postNBR);
// Scale
var dNBR = dNBR_unscaled.multiply(1000);

//==
//Layers

var maps = [];
// Map 1 Boundary + Greyscale NBR (for burnt pixel detection !!!SHADOWS = FALSE POSITIVE) + NBR + NDVI

var maptopleft = ui.Map();
maptopleft.add(ui.Label('dNBR'));

var grey = ['white', 'black'];
maptopleft.addLayer(dNBR, {min: -1000, max: 1000, palette: grey}, 'dNBR greyscale');

var sld_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap type="intervals" extended="false" >' +
      '<ColorMapEntry color="#ffffff" quantity="-500" label="-500"/>' + 
      '<ColorMapEntry color="#7a8737" quantity="-250" label="-250" />' +  
      '<ColorMapEntry color="#acbe4d" quantity="-100" label="-100" />' +
      '<ColorMapEntry color="#00ff00" quantity="100" label="100" />' +
      '<ColorMapEntry color="#fff70b" quantity="270" label="270" />' +
      '<ColorMapEntry color="#ffaf38" quantity="440" label="440" />' +
      '<ColorMapEntry color="#ff641b" quantity="660" label="660" />' +
      '<ColorMapEntry color="#a41fd6" quantity="2000" label="2000" />' +
    '</ColorMap>' +
  '</RasterSymbolizer>';

maptopleft.addLayer(dNBR.sldStyle(sld_intervals), {}, 'dNBR classified');

var thresholdlim = ee.Image([-1000, -251, -101, 99, 269, 439, 659, 2000]);
// var thresholdlim = ee.Image([-1000, -251, -101, 99, 269, 439, 659, 2000]);
var classified = dNBR.lt(thresholdlim).reduce('sum').toInt();

maptopleft.setControlVisibility(false);
maptopleft.setControlVisibility('Layers', true);
maptopleft.setOptions("SATELLITE")
maps.push(maptopleft);
maptopleft.setCenter(28.771110623071657,36.8506091081693, 7);


var visParams = {bands: ['B8', 'B4', 'B3'], max: 3048, gamma: 1};
var visParams_ndvi = {min: -0.2, max: 0.8, palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,' +
    '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'};

var image_ndvi = pre_cm_mos.normalizedDifference(['B8','B4']);

maptopleft.addLayer(image_ndvi,{min: -0.2, max: 0.8, palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,' +
    '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'}, 'NDVI' )

//--
//True color

if (platform == 'S2' | platform == 's2') {
  var vis = {bands: ['B4', 'B3', 'B2'], max: 2000, gamma: 0.70}; //lower gamma for S2
} else {
  var vis = {bands: ['B4', 'B3', 'B2'], min: 0, max: 4000, gamma: 0.70}; //^what she said
}

var left = ui.Map();
var right = ui.Map();
ui.root.clear();
ui.root.add(left);
ui.root.add(right);
ui.root.add(maptopleft);

ui.Map.Linker([left, right], 'change-bounds');

left.addLayer(pre_cm_mos, vis,'Pre-fire True Color Image - Clouds masked'); //from prev. section
right.addLayer(post_cm_mos, vis,'Post-fire True Color Image - Clouds masked');

var linker = ui.Map.Linker([ui.root.widgets().get(0), right]);
var splitPanel = ui.SplitPanel({
  firstPanel: linker.get(0),
  secondPanel: linker.get(1),
  orientation: 'horizontal',
  wipe: true,
  style: {stretch: 'both'}
});

left.add(ui.Label('Pre- & Post-fire'));
left.setControlVisibility(false);
right.add(ui.Label('Pre- & Post-fire'));
right.setControlVisibility(false);

left.setCenter(28.771110623071657,36.8506091081693, 7);
right.setCenter(28.771110623071657,36.8506091081693, 7);
left.setOptions("SATELLITE")
right.setOptions("SATELLITE")

ui.root.widgets().reset([maptopleft, splitPanel]);
ui.root.setLayout(ui.Panel.Layout.Flow('horizontal'));

//===
// Burned area
var allpix =  classified.updateMask(classified);  // mask the entire layer
var pixstats = allpix.reduceRegion({
  reducer: ee.Reducer.count(),               // count pixels in a single class
  geometry: area,
  scale: 65   //original=30 NEED TO NORMALIZE PIXEL COUNT FOR HECT + PERC
  });
var allpixels = ee.Number(pixstats.get('sum')); // extract pixel count

var arealist = [];
var areacount = function(cnr, name) {
 var singleMask =  classified.updateMask(classified.eq(cnr));  // single class
 var stats = singleMask.reduceRegion({
  reducer: ee.Reducer.count(),               // count pixels in a single class
  geometry: area,
  scale: 65     //original L8=30 NEED TO NORMALIZE PIXEL COUNT FOR HECT + PERC
  });
var pix =  ee.Number(stats.get('sum')).multiply(2.1666666667).round();
var hect = pix.multiply(900).divide(10000).multiply(2.1666666667);                // Landsat pix = 30m x 30m = 900 sqm
var perc = pix.divide(allpixels).multiply(10000).divide(100).multiply(2.1666666667);   // area percent by class rounded
arealist.push({Class: name, Pixels: pix, Hectares: hect, Percentage: perc});
};

// Severity classes 
var names2 = ['NA', 'High Severity', 'Moderate-high Severity',
'Moderate-low Severity', 'Low Severity','Unburned', 'Enhanced Regrowth, High','Enhanced Regrowth, Low'];

for (var i = 1; i < 5; i++) {
  areacount(i, names2[i]);
  }

print('Burned Area by Severity Class', arealist);

//==
// Legend

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }});

var legendTitle = ui.Label({
  value: 'Normalized Burn Ratio Classes',
  style: {fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }});
legend.add(legendTitle);
 

var makeRow = function(color, name) {
 
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,

          padding: '8px',
          margin: '0 0 4px 0'
        }});
 
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      })};
 
var palette =['7a8737', 'acbe4d', '0ae042', 'fff70b', 'ffaf38', 'ff641b', 'a41fd6', 'ffffff'];
 
var names = ['Enhanced Regrowth, High','Enhanced Regrowth, Low','Unburned', 'Low Severity',
'Moderate-low Severity', 'Moderate-high Severity', 'High Severity', 'NA'];
 
for (var i = 3; i < 7; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  

maptopleft.add(legend);
