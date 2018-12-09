
let scenewrapper;

const base_folder = './datasets/non-rigid-world/'

function start() {
  scenewrapper = new SceneWrapper();
  load_object(base_folder + 'cat0.obj', new THREE.Vector3(-100, 0, 0));
  load_object(base_folder + 'cat1.obj', new THREE.Vector3(100, 0, 0));
  display_sample('corresp-' + 'cat0', new THREE.Vector3(-100, 0, 0));
  display_sample('corresp-' + 'cat1', new THREE.Vector3(100, 0, 0));
}

function load_object(obj_file, pos) {
  const loader = new THREE.OBJLoader();
  loader.load(obj_file, (obj) => { 
    if(pos != null)
      obj.position.set(pos.x, pos.y, pos.z);
    scenewrapper.scene.add(obj); 
  });
}

function display_sample(sample_name, pos) {
  const loader = new THREE.FileLoader();
  loader.load(sample_name, (data) => {
    const lines = data.split('\n');
    console.log(lines);
    for(let line of lines) {
      const nums = parse_nums(line);
      nums[0] += pos.x;
      nums[1] += pos.y;
      nums[2] += pos.z;
      add_sphere(nums);
    }
  });
}

function parse_nums(line) {
  console.log(line);
  const str_nums = line.split(' ');
  console.log(str_nums);
  return str_nums.map(x => parseFloat(x));
}

function add_sphere(nums) {
  const [x, y, z] = nums;
  console.log(nums);
  const radius = 1;
  const resolution = 30;
  const geometry = new THREE.SphereGeometry(radius);
  const material = new THREE.MeshBasicMaterial({color : 0xEEEE00});

  const sphere = new THREE.Mesh(geometry, material);
  sphere.position.set(x, y, z);

  scenewrapper.scene.add(sphere);
}

function random_hex_color() {
  return "#000000".replace(/0/g,function(){return (~~(Math.random()*16)).toString(16);});
}

window.onload = function () {
  start();
  // display_sample();
}
