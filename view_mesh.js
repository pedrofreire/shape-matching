
let scenewrapper;

// const base_folder = 'outputs/all_runs/40';
const base_folder = 'outputs/last';

function start() {
  scenewrapper = new SceneWrapper();
  scenewrapper.renderer.setClearColor( 0xFFFFFF, 1);
  const num_colors = 200;
  const colors = [...Array(num_colors)].map(_ => random_color());

  const half_dist = 60.0;

  load_object(base_folder + '/mesh0.obj', new THREE.Vector3(-half_dist, 0, 0));
  load_object(base_folder + '/mesh1.obj', new THREE.Vector3(half_dist, 0, 0));
  display_pointset(base_folder + '/correspondence0', new THREE.Vector3(-half_dist, 0, 0), 3.0, colors);
  display_pointset(base_folder + '/correspondence1', new THREE.Vector3( half_dist, 0, 0),  3.0, colors);
  display_pointset(base_folder + '/sample0', new THREE.Vector3(-half_dist, 0, 0), 2.0);
  display_pointset(base_folder + '/sample1', new THREE.Vector3(half_dist, 0, 0), 2.0);
}

function random_color() {
  return Math.floor(Math.random() * 0xFFFFFF);
}

function load_object(obj_file, pos) {
  const loader = new THREE.OBJLoader();
  loader.load(obj_file, (obj) => { 
    if(pos != null)
      obj.position.set(pos.x, pos.y, pos.z);
    obj.children[0].material.color.setHex(0xbbbbbb);
    scenewrapper.scene.add(obj); 
  });
}

function display_pointset(sample_name, pos, radius, colors) {

  const loader = new THREE.FileLoader();
  loader.load(sample_name, (data) => {
    const lines = data.split('\n');
    for(let i = 0; i != lines.length; ++i) {
      const nums = parse_nums(lines[i]);
      nums[0] += pos.x;
      nums[1] += pos.y;
      nums[2] += pos.z;

      const color = colors != null ? colors[i] : 0xEEEE00;
      add_sphere(nums, radius, color);
    }
  });
}

function parse_nums(line) {
  const str_nums = line.split(' ');
  return str_nums.map(x => parseFloat(x));
}

function add_sphere(nums, radius, color) {
  const [x, y, z] = nums;
  const resolution = 30;
  const geometry = new THREE.SphereGeometry(radius);
  const material = new THREE.MeshBasicMaterial({color : color});

  const sphere = new THREE.Mesh(geometry, material);
  sphere.position.set(x, y, z);

  scenewrapper.scene.add(sphere);
}

window.onload = function () {
  start();
}
