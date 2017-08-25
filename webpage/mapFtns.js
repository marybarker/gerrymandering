function myMap(){
  /* CREATE MAP */
  map = new google.maps.Map(document.getElementById("mapCanvas"), {
    zoom:7, 
    center: new google.maps.LatLng(centroid[0],centroid[1]) 
  });
  
  /* ADD VTDS TO MAP */
  map.data.loadGeoJson(jsonfilename, {}, function(features){yet_to_assign = features.length;});
  
  /* SET DISTRICT OF EACH VTD TO UNASSIGNED */
  map.data.forEach(function(feature){
    feature.setProperty('District', 'Not Assigned');
  });
  
  /* COLOR VTDS */
  map.data.setStyle(function(feature){
    myColor = (feature.getProperty("District") === undefined) ? allColors['Not Assigned'] : allColors[feature.getProperty('District')];
    return ({
      fillColor: myColor, 
      strokeWeight: .1,
      fillOpacity:0.8
    });
  });
  
  /* ADD POLYGON SELECT OPTION */
  var drawingManager = new google.maps.drawing.DrawingManager({
    drawingControl: true, 
    drawingControlOptions: {
      position: google.maps.ControlPosition.TOP_CENTER, 
      drawingModes: ['polygon']
    }
  });
  drawingManager.setMap(map);
  // when polygon has been completed, add vtds inside to currentDistrict 
  google.maps.event.addListener(drawingManager, 'overlaycomplete', function(event){
    // shape of the polygon
    var theShape = event.overlay;
    // select which vtds are inside the polygon
    map.data.forEach(function(feature){
      if (feature.getProperty("District") != currentDistrict) {
        var bounds = new google.maps.LatLngBounds();
        var theGeom = feature.getGeometry();
        theGeom.forEachLatLng(function(LatLng){
          bounds.extend(LatLng);
        });
        if (google.maps.geometry.poly.containsLocation(bounds.getCenter(), theShape) ){
          if (currentDistrict != "Not Assigned") {
            yet_to_assign--;
          }
          if (feature.getProperty("District") >= 0) {
            yet_to_assign++;
          }
          feature.setProperty("District", currentDistrict);
        }
      }
    });
    // get rid of the polygon once vtds are chosen
    theShape.setMap(null);
    updateCurrentStateInfo();
  });
  
  /* CLICK SELECT OPTION */
  map.data.addListener('click', function(event){
    // get district this feature was in before click
    var value = event.feature.getProperty("District");
    // set the feature to be in district currentDistrict
    event.feature.setProperty("District", currentDistrict);
    // update total number of unassigned features.
    if (currentDistrict != "Not Assigned") {
      yet_to_assign--;
    }
    if (value >= 0) {
      yet_to_assign++;
    }
    // now update stats on current district
    updateCurrentStateInfo();
  });
  
  /* HOVER OVER VTD HIGHLIGHT OPTION */
  map.data.addListener('mouseover', function(event){
    map.data.revertStyle();
    map.data.overrideStyle(event.feature, {strokeWeight:2});
  });
  map.data.addListener('mouseout', function(event){
    map.data.revertStyle();
  });
}

