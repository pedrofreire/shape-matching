
const express = require('express');
const ejs = require('ejs');
const path = require('path');

const app = express();

app.set('view engine', 'ejs');
app.set('views', path.join(__dirname + '/'));
app.use(express.static(path.join(__dirname + '/')));

app.get('/', (req, res) => {
  const run_id = req.query.run_id ? req.query.run_id : '';
  res.render('index.ejs', {RUN_ID : run_id});
});

const port = 8088;

app.listen(port);
