package hacone.AgIn.resources

object Resources {
  // lda coefficient
  val p4c2 = List(
   +0.095469191, -0.310231519, +0.183874921, +0.136647155, -0.525053543,
   +0.234870353, +0.038720183, -0.245769328, -0.594874060, +0.349912265, -0.739901459,
   +0.017933386, -0.212674876, -0.037288583, -0.002777067, -0.046708805,
   +0.011699703, +0.016603473, +0.048500211, +0.003669691, +0.022911584).map(-_)

  val p5c3 = List(
  -0.080311809, 0.1541256439, -0.1726490191, -0.0753845871, 0.4206795098,
  -0.1237638162, -0.0514447233, 0.1744636532, 0.4224327126, -0.5720887972,
  0.3828247197,
  0.0174045372, 0.2223319501, -0.0125660944, -0.0171837688, -0.0080998961,
  -0.0501252467, 0.0224256831, -0.0493460137, 0.0296150985, 0.0078371746)
  
}
