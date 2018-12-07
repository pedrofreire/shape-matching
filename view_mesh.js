

function start() {
  const scenewrapper = new SceneWrapper();
  const loader = new THREE.OBJLoader();
  loader.load('./datasets/non-rigid-world/victoria0.obj', (obj) => { 
  // loader.load('./datasets/simple/reg_tetra.obj', (obj) => { 
    scenewrapper.scene.add(obj); 
  });
}

window.onload = function () {
  start();
}
